/// Gauss Lemniscate (ϖ) ve Eliptik Sabitler İçin Mod-4 Hipergeometrik PCF Tanımları
/// Kuyruk Kesir (Tail Evaluation: n=1..N) Mimarisi
/// Intel Core i3-N305 (8 Gracemont E-Core, AVX2, 6MB L3) için optimize edilmiştir.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PcfFamily {
    /// Aile 1: Asimetrik Çarpanlı Lineer Gauss-Pfaff Ailesi
    /// a_n = c1*n + c2
    /// b_n = -c3 * (p1*n + q1) * (p2*n + q2)
    LinearGaussPfaff,

    /// Aile 2: Kuadratik a_n ve 4. Derece Mod-4 Çarpanlı b_n (Genelleştirilmiş Contiguous)
    /// a_n = c1*n^2 + c2*n + c3
    /// b_n = -c4 * n * (n + k) * (p1*n + q1) * (p2*n + q2)
    QuadraticContiguous,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PcfCandidate {
    pub family: PcfFamily,
    // a_n katsayıları
    pub c1: i64,
    pub c2: i64,
    pub c3: i64, // Aile 1'de 0, Aile 2'de c3
    // b_n katsayıları
    pub c_factor: i64, // Aile 1'de c3, Aile 2'de c4
    pub k: i64,        // Aile 1'de 0, Aile 2'de k in [0, 2]
    pub p1: i64,
    pub q1: i64,
    pub p2: i64,
    pub q2: i64,
}

pub fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

impl PcfCandidate {
    /// a_n matematiksel string gösterimi
    pub fn a_string(&self) -> String {
        match self.family {
            PcfFamily::LinearGaussPfaff => {
                let base = if self.c1 == 1 {
                    "n".to_string()
                } else {
                    format!("{}n", self.c1)
                };

                if self.c2 > 0 {
                    format!("{} + {}", base, self.c2)
                } else if self.c2 < 0 {
                    format!("{} - {}", base, -self.c2)
                } else {
                    base
                }
            }
            PcfFamily::QuadraticContiguous => {
                let mut parts = Vec::with_capacity(3);
                if self.c1 == 1 {
                    parts.push("n^2".to_string());
                } else {
                    parts.push(format!("{}*n^2", self.c1));
                }

                if self.c2 > 0 {
                    if self.c2 == 1 {
                        parts.push("n".to_string());
                    } else {
                        parts.push(format!("{}*n", self.c2));
                    }
                } else if self.c2 < 0 {
                    if self.c2 == -1 {
                        parts.push("- n".to_string());
                    } else {
                        parts.push(format!("- {}*n", -self.c2));
                    }
                }

                if self.c3 > 0 {
                    parts.push(format!("{}", self.c3));
                } else if self.c3 < 0 {
                    parts.push(format!("- {}", -self.c3));
                }

                parts.join(" + ").replace("+ - ", "- ")
            }
        }
    }

    /// b_n matematiksel string gösterimi
    pub fn b_string(&self) -> String {
        let fmt_factor = |p: i64, q: i64| -> String {
            let p_str = if p == 1 {
                "n".to_string()
            } else {
                format!("{}n", p)
            };
            if q > 0 {
                format!("({} + {})", p_str, q)
            } else if q < 0 {
                format!("({} - {})", p_str, -q)
            } else {
                p_str
            }
        };

        let f1_str = fmt_factor(self.p1, self.q1);
        let f2_str = fmt_factor(self.p2, self.q2);

        match self.family {
            PcfFamily::LinearGaussPfaff => {
                let prefix = if self.c_factor == 1 {
                    "-".to_string()
                } else {
                    format!("-{} * ", self.c_factor)
                };

                // Eğer f1_str parantezsiz "n" ise ve prefix "-" ise "-n * ..."
                format!("{}{} * {}", prefix, f1_str, f2_str)
            }
            PcfFamily::QuadraticContiguous => {
                let prefix = if self.c_factor == 1 {
                    "-".to_string()
                } else {
                    format!("-{} * ", self.c_factor)
                };

                let k_str = if self.k == 0 {
                    "n^2".to_string()
                } else {
                    format!("n * (n + {})", self.k)
                };

                format!("{}{} * {} * {}", prefix, k_str, f1_str, f2_str)
            }
        }
    }

    /// Aynı değere yakınsayan skaler katları (a_n -> d*a_n, b_n -> d^2*b_n) tespit eder
    #[inline(always)]
    pub fn is_redundant_scalar_multiple(&self) -> bool {
        match self.family {
            PcfFamily::LinearGaussPfaff => {
                let g_a = gcd(self.c1, self.c2.abs());
                if g_a <= 1 {
                    return false;
                }
                // b_n katsayıları: c_factor * (p1*p2*n^2 + (p1*q2 + p2*q1)*n + q1*q2)
                let a_coef = self.c_factor * self.p1 * self.p2;
                let b_coef = self.c_factor * (self.p1 * self.q2 + self.p2 * self.q1);
                let c_coef = self.c_factor * self.q1 * self.q2;
                let g_b = gcd(a_coef, gcd(b_coef.abs(), c_coef.abs()));

                // g_a'nın herhangi bir d > 1 böleni için d^2 | g_b ise skaler kattır
                let mut d = 2;
                while d * d <= g_a {
                    if g_a % d == 0 {
                        if g_b % (d * d) == 0 {
                            return true;
                        }
                        let other = g_a / d;
                        if other > 1 && g_b % (other * other) == 0 {
                            return true;
                        }
                    }
                    d += 1;
                }
                if g_a > 1 && g_b % (g_a * g_a) == 0 {
                    return true;
                }
                false
            }
            PcfFamily::QuadraticContiguous => {
                let g_a = gcd(self.c1, gcd(self.c2.abs(), self.c3.abs()));
                if g_a <= 1 {
                    return false;
                }
                // Basit kontrol: c_factor g_a^2 ile tam bölünebiliyorsa
                if self.c_factor % (g_a * g_a) == 0 {
                    return true;
                }
                false
            }
        }
    }
}

/// Aşama 1: O(1) Matematiksel Ön-Eleme (Pre-filter)
#[inline(always)]
pub fn prefilter(candidate: &PcfCandidate) -> bool {
    match candidate.family {
        PcfFamily::LinearGaussPfaff => {
            if candidate.c1 <= 0 || candidate.c_factor <= 0 {
                return false;
            }
            // Dejenere b_1 == 0 eleme
            if candidate.p1 + candidate.q1 == 0 || candidate.p2 + candidate.q2 == 0 {
                return false;
            }
            // Darboux limiti: 4 * c3 * p1 * p2 <= c1^2
            let c1_sq = candidate.c1 * candidate.c1;
            if 4 * candidate.c_factor * candidate.p1 * candidate.p2 > c1_sq {
                return false;
            }
            // a_1 = c1 + c2 kontrolü
            if candidate.c1 + candidate.c2 <= 0 {
                return false;
            }
            // Skaler kat indirgemesi
            if candidate.is_redundant_scalar_multiple() {
                return false;
            }
            true
        }
        PcfFamily::QuadraticContiguous => {
            if candidate.c1 <= 0 || candidate.c_factor <= 0 {
                return false;
            }
            // Dejenere b_1 == 0 eleme
            if candidate.p1 + candidate.q1 == 0 || candidate.p2 + candidate.q2 == 0 {
                return false;
            }
            // Darboux limiti: 4 * c4 * p1 * p2 <= c1^2
            let c1_sq = candidate.c1 * candidate.c1;
            if 4 * candidate.c_factor * candidate.p1 * candidate.p2 > c1_sq {
                return false;
            }
            // a_1 = c1 + c2 + c3 kontrolü
            if candidate.c1 + candidate.c2 + candidate.c3 <= 0 {
                return false;
            }
            if candidate.is_redundant_scalar_multiple() {
                return false;
            }
            true
        }
    }
}

/// Aşama 2: f64 Wallis Matris İterasyonu (Fast Float Filter)
/// KUYRUK KESİR (Tail Evaluation): a_0 terimi OLMADAN doğrudan n=1..35 hesaplar.
///
/// Euler-Wallis rekürans başlangıcı:
/// A_0 = 0.0, B_0 = 1.0
/// A_1 = b_1, B_1 = a_1
/// n = 2..35 için:
///   A_n = a_n * A_{n-1} + b_n * A_{n-2}
///   B_n = a_n * B_{n-1} + b_n * B_{n-2}
#[inline(always)]
pub fn evaluate_wallis_35(candidate: &PcfCandidate) -> Option<f64> {
    // 1. Adım katsayıları
    let (a1, b1) = match candidate.family {
        PcfFamily::LinearGaussPfaff => {
            let a1 = (candidate.c1 + candidate.c2) as f64;
            let b1 = (-candidate.c_factor * (candidate.p1 + candidate.q1) * (candidate.p2 + candidate.q2)) as f64;
            (a1, b1)
        }
        PcfFamily::QuadraticContiguous => {
            let a1 = (candidate.c1 + candidate.c2 + candidate.c3) as f64;
            let b1 = (-candidate.c_factor * (1 + candidate.k) * (candidate.p1 + candidate.q1) * (candidate.p2 + candidate.q2)) as f64;
            (a1, b1)
        }
    };

    if b1 == 0.0 || a1 == 0.0 {
        return None;
    }

    let mut a_prev2 = 0.0f64; // A_0
    let mut a_prev1 = b1;     // A_1
    let mut b_prev2 = 1.0f64; // B_0
    let mut b_prev1 = a1;     // B_1

    let c1 = candidate.c1 as f64;
    let c2 = candidate.c2 as f64;
    let c3 = candidate.c3 as f64;
    let c_factor = candidate.c_factor as f64;
    let k = candidate.k as f64;
    let p1 = candidate.p1 as f64;
    let q1 = candidate.q1 as f64;
    let p2 = candidate.p2 as f64;
    let q2 = candidate.q2 as f64;

    for n_int in 2..=35i32 {
        let n = n_int as f64;

        let (an, bn) = match candidate.family {
            PcfFamily::LinearGaussPfaff => {
                let an = c1 * n + c2;
                let bn = -c_factor * (p1 * n + q1) * (p2 * n + q2);
                (an, bn)
            }
            PcfFamily::QuadraticContiguous => {
                let an = (c1 * n + c2) * n + c3;
                let bn = -c_factor * n * (n + k) * (p1 * n + q1) * (p2 * n + q2);
                (an, bn)
            }
        };

        let a_curr = an * a_prev1 + bn * a_prev2;
        let b_curr = an * b_prev1 + bn * b_prev2;

        if !a_curr.is_finite() || !b_curr.is_finite() {
            return None;
        }

        a_prev2 = a_prev1;
        a_prev1 = a_curr;
        b_prev2 = b_prev1;
        b_prev1 = b_curr;

        // Taşma ve yetersizlik önleyici ölçekleme
        if b_curr.abs() > 1.0e80 {
            a_prev1 *= 1.0e-80;
            a_prev2 *= 1.0e-80;
            b_prev1 *= 1.0e-80;
            b_prev2 *= 1.0e-80;
        } else if n_int % 10 == 0 && b_curr.abs() > 1.0e40 {
            a_prev1 *= 1.0e-40;
            a_prev2 *= 1.0e-40;
            b_prev1 *= 1.0e-40;
            b_prev2 *= 1.0e-40;
        }
    }

    if b_prev1.abs() < 1.0e-300 {
        return None;
    }

    let val = a_prev1 / b_prev1;
    if val.is_finite() {
        Some(val)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test 1: a_n = 4n + 1, b_n = -n(2n - 1) => K = 2/ϖ - 1 ≈ -0.2372402364981868
    #[test]
    fn test_known_lemniscate_pcf_1() {
        let cand = PcfCandidate {
            family: PcfFamily::LinearGaussPfaff,
            c1: 4,
            c2: 1,
            c3: 0,
            c_factor: 1,
            k: 0,
            p1: 1,
            q1: 0,
            p2: 2,
            q2: -1,
        };

        assert!(prefilter(&cand), "Test 1 prefilter başarısız oldu!");

        let res = evaluate_wallis_35(&cand);
        assert!(res.is_some(), "Test 1 evaluate_wallis_35 None döndü!");

        let val = res.unwrap();
        // ϖ ≈ 2.6220575542921198
        let expected = (2.0 / 2.6220575542921198) - 1.0; // ≈ -0.2372402364981868
        let diff = (val - expected).abs();
        assert!(
            diff < 1.0e-12,
            "Test 1 değeri beklenen ile uyuşmuyor! Hesaplanan: {:.16}, Beklenen: {:.16}, Fark: {:e}",
            val,
            expected,
            diff
        );
    }

    /// Test 2: a_n = 6n, b_n = -(2n - 1)(4n - 3) => K = 2/ϖ - 1 ≈ -0.2372402364981868
    #[test]
    fn test_known_lemniscate_pcf_2() {
        let cand = PcfCandidate {
            family: PcfFamily::LinearGaussPfaff,
            c1: 6,
            c2: 0,
            c3: 0,
            c_factor: 1,
            k: 0,
            p1: 2,
            q1: -1,
            p2: 4,
            q2: -3,
        };

        assert!(prefilter(&cand), "Test 2 prefilter başarısız oldu!");

        let res = evaluate_wallis_35(&cand);
        assert!(res.is_some(), "Test 2 evaluate_wallis_35 None döndü!");

        let val = res.unwrap();
        let expected = (2.0 / 2.6220575542921198) - 1.0; // ≈ -0.2372402364981868
        let diff = (val - expected).abs();
        assert!(
            diff < 1.0e-9,
            "Test 2 değeri beklenen ile uyuşmuyor! Hesaplanan: {:.16}, Beklenen: {:.16}, Fark: {:e}",
            val,
            expected,
            diff
        );
    }
}
