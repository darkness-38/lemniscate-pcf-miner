use crate::pcf::{PcfCandidate, PcfFamily};
use dashu::base::Abs;
use dashu::float::round::mode::HalfAway;
use dashu::float::{Context, DBig};
use dashu::integer::IBig;
use std::fs::OpenOptions;
use std::io::Write;
use std::str::FromStr;

/// 1. Gauss Lemniscate Sabiti (ϖ)
pub const LEMNISCATE_STR: &str =
    "2.62205755429211981046483958989111941368275495143162316281682170380079058707041425406540817080011375";
pub const LEMNISCATE_F64: f64 = 2.6220575542921198;

/// 2. Gauss AGM Sabiti (pi / ϖ)
pub const GAUSS_AGM_STR: &str =
    "1.1981402347355922074399224922803238782272126632156533";
pub const GAUSS_AGM_F64: f64 = 1.1981402347355922;

/// 3. İkinci Lemniscate / Eliptik Sabit K(1/√2) (ϖ / √2)
pub const ELLIPTIC_K_STR: &str =
    "1.8540746773013719184338503471952600462175988235219845";
pub const ELLIPTIC_K_F64: f64 = 1.854074677301372;

/// 4. Karesel Lemniscate Oranı (ϖ² / pi)
pub const QUADRATIC_RATIO_STR: &str =
    "2.1884396152437628287383794354746811231802927236592264";
pub const QUADRATIC_RATIO_F64: f64 = 2.188439615243763;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BaseConstant {
    Lemniscate,      // ϖ
    GaussAgm,        // pi / ϖ
    EllipticK,       // ϖ / sqrt(2)
    QuadraticRatio,  // ϖ^2 / pi
}

#[derive(Clone, Debug)]
pub struct TargetEntry {
    pub value_f64: f64,
    pub name: String,
    #[allow(dead_code)]
    pub formula: String,
    pub base: BaseConstant,
    /// Möbius katsayıları: Target = (m1*S + m2) / (m3*S + m4)
    pub m: (i64, i64, i64, i64),
}

#[derive(Clone, Debug)]
pub struct VerificationResult {
    pub target_name: String,
    pub a_repr: String,
    pub b_repr: String,
    pub value_str: String,
    pub error_str: String,
    pub steps: usize,
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

fn gcd4(a: i64, b: i64, c: i64, d: i64) -> i64 {
    gcd(gcd(a, b), gcd(c, d))
}

fn add_targets_for_constant(
    targets: &mut Vec<TargetEntry>,
    base: BaseConstant,
    val_f64: f64,
    symbol: &str,
    full_name: &str,
    depth: i64,
) {
    // 1. Doğrudan sabit
    targets.push(TargetEntry {
        value_f64: val_f64,
        name: format!("{} ({})", symbol, full_name),
        formula: symbol.to_string(),
        base,
        m: (1, 0, 0, 1),
    });

    // 2. Möbius Dönüşümleri: (m1*S + m2) / (m3*S + m4)
    for m1 in -depth..=depth {
        for m2 in -depth..=depth {
            for m3 in -depth..=depth {
                for m4 in -depth..=depth {
                    let det = m1 * m4 - m2 * m3;
                    if det == 0 {
                        continue;
                    }

                    let (s1, s2, s3, s4) = if m3 > 0 || (m3 == 0 && m4 > 0) {
                        (m1, m2, m3, m4)
                    } else {
                        (-m1, -m2, -m3, -m4)
                    };

                    let g = gcd4(s1.abs(), s2.abs(), s3.abs(), s4.abs());
                    if g > 1 {
                        continue;
                    }

                    let denom_f64 = (s3 as f64) * val_f64 + (s4 as f64);
                    if denom_f64.abs() < 1.0e-4 {
                        continue;
                    }

                    let num_f64 = (s1 as f64) * val_f64 + (s2 as f64);
                    let target_val = num_f64 / denom_f64;

                    if !target_val.is_finite() {
                        continue;
                    }

                    if s1 == 1 && s2 == 0 && s3 == 0 && s4 == 1 {
                        continue;
                    }

                    let (name, formula) = format_mobius_name(s1, s2, s3, s4, symbol, full_name);

                    targets.push(TargetEntry {
                        value_f64: target_val,
                        name,
                        formula,
                        base,
                        m: (s1, s2, s3, s4),
                    });
                }
            }
        }
    }
}

fn format_mobius_name(
    m1: i64,
    m2: i64,
    m3: i64,
    m4: i64,
    sym: &str,
    full_name: &str,
) -> (String, String) {
    if m1 == 0 && m2 == 1 && m3 == 1 && m4 == 0 {
        return (
            format!("1/{} (Reciprocal {})", sym, full_name),
            format!("1/{}", sym),
        );
    }
    if m3 == 0 && m4 == 1 {
        let f = if m2 == 0 {
            format!("{}*{}", m1, sym)
        } else if m2 > 0 {
            format!("{}*{} + {}", m1, sym, m2)
        } else {
            format!("{}*{} - {}", m1, sym, -m2)
        };
        return (format!("{} [Linear {}]", f, sym), f);
    }
    if m1 == 0 {
        let denom = if m4 == 0 {
            format!("{}*{}", m3, sym)
        } else if m4 > 0 {
            format!("{}*{} + {}", m3, sym, m4)
        } else {
            format!("{}*{} - {}", m3, sym, -m4)
        };
        let f = format!("{}/({})", m2, denom);
        return (format!("{} [Möbius Inverse {}]", f, sym), f);
    }

    let num = if m2 == 0 {
        format!("{}*{}", m1, sym)
    } else if m2 > 0 {
        format!("{}*{} + {}", m1, sym, m2)
    } else {
        format!("{}*{} - {}", m1, sym, -m2)
    };

    let denom = if m4 == 0 {
        format!("{}*{}", m3, sym)
    } else if m4 > 0 {
        format!("{}*{} + {}", m3, sym, m4)
    } else {
        format!("{}*{} - {}", m3, sym, -m4)
    };

    let f = format!("({}) / ({})", num, denom);
    (format!("{} [{} Möbius]", f, sym), f)
}

/// 4 temel sabit ve [-depth, depth] Möbius dönüşümlerinden hedef tabloyu oluşturur.
pub fn build_target_table(depth: i64) -> Vec<TargetEntry> {
    let mut targets = Vec::new();

    // 1. Lemniscate Sabiti (ϖ)
    add_targets_for_constant(
        &mut targets,
        BaseConstant::Lemniscate,
        LEMNISCATE_F64,
        "ϖ",
        "Lemniscate",
        depth,
    );

    // 2. Gauss AGM Sabiti (pi / ϖ)
    add_targets_for_constant(
        &mut targets,
        BaseConstant::GaussAgm,
        GAUSS_AGM_F64,
        "(pi/ϖ)",
        "Gauss AGM",
        depth,
    );

    // 3. Tam Eliptik İntegral K(1/√2) (ϖ / √2)
    add_targets_for_constant(
        &mut targets,
        BaseConstant::EllipticK,
        ELLIPTIC_K_F64,
        "K(1/√2)",
        "Elliptic K",
        depth,
    );

    // 4. Karesel Oran (ϖ² / pi)
    add_targets_for_constant(
        &mut targets,
        BaseConstant::QuadraticRatio,
        QUADRATIC_RATIO_F64,
        "(ϖ²/pi)",
        "Quadratic Ratio",
        depth,
    );

    // İkili arama için sırala ve tekilleştir
    targets.sort_by(|a, b| a.value_f64.partial_cmp(&b.value_f64).unwrap());
    targets.dedup_by(|a, b| (a.value_f64 - b.value_f64).abs() < 1.0e-12);
    targets
}

/// Sıralı hedef tablosunda val değerine tolerans dahilinde eşleşen hedefleri bulur.
#[inline(always)]
pub fn find_matching_targets<'a>(
    val: f64,
    targets: &'a [TargetEntry],
    tolerance: f64,
) -> Vec<&'a TargetEntry> {
    let mut matches = Vec::new();

    let low_idx = targets.partition_point(|t| t.value_f64 < val - tolerance);
    for target in &targets[low_idx..] {
        if target.value_f64 > val + tolerance {
            break;
        }
        if (target.value_f64 - val).abs() < tolerance {
            matches.push(target);
        }
    }

    matches
}

/// Aşama 3: Yüksek Duyarlıklı Doğrulama (Arbitrary-Precision Verification)
/// KUYRUK KESİR (Tail Evaluation): n=1..150 adım tam sayı (IBig) reküransı
/// A_0 = 0, B_0 = 1
/// A_1 = b_1, B_1 = a_1
/// n=2..150 için:
///   A_n = a_n * A_{n-1} + b_n * A_{n-2}
///   B_n = a_n * B_{n-1} + b_n * B_{n-2}
/// En az 40 basamak tam sayı (dashu) eşleşmesi arar (|V_150 - Target| < 10^-40).
pub fn verify_stage3(candidate: &PcfCandidate, target: &TargetEntry) -> Option<VerificationResult> {
    const STEPS: usize = 150;
    const PRECISION_DIGITS: usize = 160;

    let (a1_val, b1_val) = match candidate.family {
        PcfFamily::LinearGaussPfaff => {
            let a1 = candidate.c1 + candidate.c2;
            let b1 = -candidate.c_factor * (candidate.p1 + candidate.q1) * (candidate.p2 + candidate.q2);
            (a1, b1)
        }
        PcfFamily::QuadraticContiguous => {
            let a1 = candidate.c1 + candidate.c2 + candidate.c3;
            let b1 = -candidate.c_factor * (1 + candidate.k) * (candidate.p1 + candidate.q1) * (candidate.p2 + candidate.q2);
            (a1, b1)
        }
    };

    if b1_val == 0 || a1_val == 0 {
        return None;
    }

    let mut a_prev2 = IBig::from(0);       // A_0
    let mut a_prev1 = IBig::from(b1_val);  // A_1
    let mut b_prev2 = IBig::from(1);       // B_0
    let mut b_prev1 = IBig::from(a1_val);  // B_1

    for n in 2..=STEPS as i64 {
        let an_val = match candidate.family {
            PcfFamily::LinearGaussPfaff => candidate.c1 * n + candidate.c2,
            PcfFamily::QuadraticContiguous => {
                candidate.c1 * n * n + candidate.c2 * n + candidate.c3
            }
        };
        let an_ibig = IBig::from(an_val);

        let bn_val = match candidate.family {
            PcfFamily::LinearGaussPfaff => {
                -candidate.c_factor * (candidate.p1 * n + candidate.q1) * (candidate.p2 * n + candidate.q2)
            }
            PcfFamily::QuadraticContiguous => {
                -candidate.c_factor * n * (n + candidate.k) * (candidate.p1 * n + candidate.q1) * (candidate.p2 * n + candidate.q2)
            }
        };
        let bn_ibig = IBig::from(bn_val);

        let a_curr = &an_ibig * &a_prev1 + &bn_ibig * &a_prev2;
        let b_curr = &an_ibig * &b_prev1 + &bn_ibig * &b_prev2;

        a_prev2 = a_prev1;
        a_prev1 = a_curr;
        b_prev2 = b_prev1;
        b_prev1 = b_curr;
    }

    if b_prev1 == IBig::from(0) {
        return None;
    }

    let ctx: Context<HalfAway> = Context::new(PRECISION_DIGITS);
    let d_num = DBig::from_parts(a_prev1, 0);
    let d_den = DBig::from_parts(b_prev1, 0);

    // V_150 = A_150 / B_150
    let v_150 = ctx.div(&d_num.repr(), &d_den.repr()).value();

    // İlgili temel sabiti seç
    let base_str = match target.base {
        BaseConstant::Lemniscate => LEMNISCATE_STR,
        BaseConstant::GaussAgm => GAUSS_AGM_STR,
        BaseConstant::EllipticK => ELLIPTIC_K_STR,
        BaseConstant::QuadraticRatio => QUADRATIC_RATIO_STR,
    };

    let base_exact = DBig::from_str(base_str).ok()?;
    let (m1, m2, m3, m4) = target.m;

    let target_exact = if m1 == 1 && m2 == 0 && m3 == 0 && m4 == 1 {
        base_exact
    } else {
        let dm1 = DBig::from(m1);
        let dm2 = DBig::from(m2);
        let dm3 = DBig::from(m3);
        let dm4 = DBig::from(m4);

        let num = ctx.mul(&dm1.repr(), &base_exact.repr()).value();
        let num = ctx.add(&num.repr(), &dm2.repr()).value();

        let den = ctx.mul(&dm3.repr(), &base_exact.repr()).value();
        let den = ctx.add(&den.repr(), &dm4.repr()).value();

        if den == DBig::from(0) {
            return None;
        }

        ctx.div(&num.repr(), &den.repr()).value()
    };

    // Fark: |V_150 - Target|
    let diff = (&v_150 - &target_exact).abs();

    // 40 basamak tolerans eşiği: 10^-40
    let threshold = DBig::from_str("0.0000000000000000000000000000000000000001").unwrap();

    if diff < threshold {
        let value_full = v_150.to_string();
        let display_val = if value_full.len() > 24 {
            format!("{}...", &value_full[..24])
        } else {
            value_full
        };

        let err_str = format_scientific(&diff);

        Some(VerificationResult {
            target_name: target.name.clone(),
            a_repr: candidate.a_string(),
            b_repr: candidate.b_string(),
            value_str: display_val,
            error_str: err_str,
            steps: STEPS,
        })
    } else {
        None
    }
}

pub fn format_scientific(d: &DBig) -> String {
    let s = d.to_string();
    if s == "0" || s == "0.0" {
        return "0.00e0".to_string();
    }

    if let Some(rest) = s.strip_prefix("0.") {
        let mut zeros = 0usize;
        let mut chars = rest.chars();
        while let Some(c) = chars.next() {
            if c == '0' {
                zeros += 1;
            } else {
                let d0 = c;
                let d1 = chars.next().unwrap_or('0');
                let exp = zeros + 1;
                return format!("{}.{}e-{}", d0, d1, exp);
            }
        }
        return "0.00e0".to_string();
    }

    s
}

pub fn save_identity(res: &VerificationResult, filename: &str) -> std::io::Result<()> {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(filename)?;

    writeln!(file, "[BULDUM!]")?;
    writeln!(file, "Hedef     : {}", res.target_name)?;
    writeln!(file, "a_n       : {}", res.a_repr)?;
    writeln!(file, "b_n       : {}", res.b_repr)?;
    writeln!(file, "Değer     : {}", res.value_str)?;
    writeln!(file, "Hata Payı : {}", res.error_str)?;
    writeln!(file, "Adım      : {} iterasyon\n", res.steps)?;
    file.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_target_table_construction() {
        let table = build_target_table(2);
        assert!(!table.is_empty());
        for i in 1..table.len() {
            assert!(
                table[i].value_f64 >= table[i - 1].value_f64,
                "Tablo sıralı değil!"
            );
        }
    }

    #[test]
    fn test_find_matching_targets() {
        let table = build_target_table(2);
        let matches = find_matching_targets(LEMNISCATE_F64, &table, 1.0e-5);
        assert!(!matches.is_empty(), "Lemniscate sabiti bulunamadı!");
    }

    #[test]
    fn test_stage3_verification_known_test1() {
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

        let target = TargetEntry {
            value_f64: (2.0 / LEMNISCATE_F64) - 1.0,
            name: "(2 - ϖ) / ϖ [ϖ Möbius]".to_string(),
            formula: "(2 - ϖ) / ϖ".to_string(),
            base: BaseConstant::Lemniscate,
            m: (-1, 2, 1, 0),
        };

        let verified = verify_stage3(&cand, &target);
        assert!(verified.is_some(), "Test 1 Stage 3 doğrulaması başarısız!");
        let res = verified.unwrap();
        println!("Test 1 Doğrulandı: {}", res.error_str);
    }

    #[test]
    fn test_stage3_verification_known_test2() {
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

        let target = TargetEntry {
            value_f64: (2.0 / LEMNISCATE_F64) - 1.0,
            name: "(2 - ϖ) / ϖ [ϖ Möbius]".to_string(),
            formula: "(2 - ϖ) / ϖ".to_string(),
            base: BaseConstant::Lemniscate,
            m: (-1, 2, 1, 0),
        };

        let verified = verify_stage3(&cand, &target);
        assert!(verified.is_some(), "Test 2 Stage 3 doğrulaması başarısız!");
        let res = verified.unwrap();
        println!("Test 2 Doğrulandı: {}", res.error_str);
    }
}
