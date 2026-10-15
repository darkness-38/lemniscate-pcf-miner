use crate::verify::{
    find_matching_targets, format_scientific, BaseConstant, TargetEntry, ELLIPTIC_K_STR,
    GAUSS_AGM_STR, LEMNISCATE_STR, QUADRATIC_RATIO_STR,
};
use dashu::base::Abs;
use dashu::float::round::mode::HalfAway;
use dashu::float::{Context, DBig};
use dashu::integer::IBig;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

#[derive(Clone, Debug)]
pub struct RawIdentity {
    pub target: String,
    pub a_str: String,
    pub b_str: String,
    pub val_str: String,
    pub err_str: String,
    pub steps: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct CanonicalIdentity {
    pub target: String,
    pub root_a: String,
    pub root_b: String,
    pub representative_val: String,
    pub best_err: String,
    pub count: usize,
}

/// discovered_lemniscate_identities.txt dosyasını okuyup ayrıştırır
pub fn parse_discovered_identities(path: &str) -> Vec<RawIdentity> {
    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Dosya açılamadı {}: {}", path, e);
            return Vec::new();
        }
    };

    let reader = BufReader::new(file);
    let mut identities = Vec::new();

    let mut current_target = String::new();
    let mut current_a = String::new();
    let mut current_b = String::new();
    let mut current_val = String::new();
    let mut current_err = String::new();
    let mut current_steps = 150;

    for line in reader.lines().flatten() {
        let line = line.trim();
        if line == "[BULDUM!]" {
            current_target.clear();
            current_a.clear();
            current_b.clear();
            current_val.clear();
            current_err.clear();
        } else if let Some(t) = line.strip_prefix("Hedef     : ") {
            current_target = t.trim().to_string();
        } else if let Some(a) = line.strip_prefix("a_n       : ") {
            current_a = a.trim().to_string();
        } else if let Some(b) = line.strip_prefix("b_n       : ") {
            current_b = b.trim().to_string();
        } else if let Some(v) = line.strip_prefix("Değer     : ") {
            current_val = v.trim().to_string();
        } else if let Some(e) = line.strip_prefix("Hata Payı : ") {
            current_err = e.trim().to_string();
        } else if let Some(s) = line.strip_prefix("Adım      : ") {
            if let Some(num) = s.split_whitespace().next() {
                current_steps = num.parse().unwrap_or(150);
            }
            if !current_target.is_empty() && !current_a.is_empty() && !current_b.is_empty() {
                identities.push(RawIdentity {
                    target: current_target.clone(),
                    a_str: current_a.clone(),
                    b_str: current_b.clone(),
                    val_str: current_val.clone(),
                    err_str: current_err.clone(),
                    steps: current_steps,
                });
            }
        }
    }

    identities
}

/// İndirgeme ve kök kimlik analizi:
/// 103 kimliği değerlerine ve çarpan sadeleşmelerine göre gruplar.
pub fn reduce_identities(raw_list: &[RawIdentity]) -> Vec<CanonicalIdentity> {
    // 1. Değerlerine göre grupla (virgülden sonra ilk 18 hane)
    let mut groups: HashMap<String, Vec<&RawIdentity>> = HashMap::new();
    for id in raw_list {
        let prefix = if id.val_str.len() >= 20 {
            &id.val_str[..20]
        } else {
            &id.val_str
        };
        groups.entry(prefix.to_string()).or_default().push(id);
    }

    let mut canonical_list = Vec::new();

    for (_val_key, ids) in groups {
        // En basit (en kısa / lineer a_n) olanı temsilci kök kimlik seç
        let mut best_id = ids[0];
        let mut min_degree = 999;
        let mut min_len = 9999;

        for id in &ids {
            let deg = if id.a_str.contains("n^2") { 2 } else { 1 };
            let len = id.a_str.len() + id.b_str.len();
            if deg < min_degree || (deg == min_degree && len < min_len) {
                min_degree = deg;
                min_len = len;
                best_id = id;
            }
        }

        // En iyi hata payını bul
        let mut best_err = best_id.err_str.clone();
        for id in &ids {
            if parse_err_exp(&id.err_str) > parse_err_exp(&best_err) {
                best_err = id.err_str.clone();
            }
        }

        canonical_list.push(CanonicalIdentity {
            target: best_id.target.clone(),
            root_a: best_id.a_str.clone(),
            root_b: best_id.b_str.clone(),
            representative_val: best_id.val_str.clone(),
            best_err,
            count: ids.len(),
        });
    }

    // Sırala: En yüksek hassasiyet (en küçük hata payı / en büyük üs) en başta
    canonical_list.sort_by(|a, b| {
        parse_err_exp(&b.best_err)
            .cmp(&parse_err_exp(&a.best_err))
            .then_with(|| a.root_a.cmp(&b.root_a))
    });

    canonical_list
}

fn parse_err_exp(err: &str) -> i32 {
    if let Some(pos) = err.find("e-") {
        err[pos + 2..].parse::<i32>().unwrap_or(0)
    } else {
        0
    }
}

/// Verilen an(n) ve bn(n) kapalı fonksiyonlarını 150 adım arbitrary precision ile hesaplar
pub fn evaluate_custom_pcf<FA, FB>(
    an_fn: FA,
    bn_fn: FB,
    steps: usize,
    precision_digits: usize,
) -> Option<DBig>
where
    FA: Fn(i64) -> i64,
    FB: Fn(i64) -> i64,
{
    let a1 = an_fn(1);
    let b1 = bn_fn(1);
    if a1 == 0 || b1 == 0 {
        return None;
    }

    let mut a_prev2 = IBig::from(0);
    let mut a_prev1 = IBig::from(b1);
    let mut b_prev2 = IBig::from(1);
    let mut b_prev1 = IBig::from(a1);

    for n in 2..=steps as i64 {
        let an = IBig::from(an_fn(n));
        let bn = IBig::from(bn_fn(n));

        let a_curr = &an * &a_prev1 + &bn * &a_prev2;
        let b_curr = &an * &b_prev1 + &bn * &b_prev2;

        a_prev2 = a_prev1;
        a_prev1 = a_curr;
        b_prev2 = b_prev1;
        b_prev1 = b_curr;
    }

    if b_prev1 == IBig::from(0) {
        return None;
    }

    let ctx: Context<HalfAway> = Context::new(precision_digits);
    let d_num = DBig::from_parts(a_prev1, 0);
    let d_den = DBig::from_parts(b_prev1, 0);
    Some(ctx.div(&d_num.repr(), &d_den.repr()).value())
}

#[derive(Clone, Debug)]
pub struct ParametricResult {
    pub k: i64,
    pub a_str: String,
    pub b_str: String,
    pub value_str: String,
    pub target_name: String,
    pub target_formula: String,
    pub error_str: String,
    pub matching_digits: usize,
}

/// Aile A ve Aile B parametrik analiz fonksiyonu
pub fn run_parametric_analysis(
    targets: &[TargetEntry],
) -> (
    Vec<ParametricResult>,
    Vec<ParametricResult>,
    Vec<ParametricResult>,
) {
    let ctx: Context<HalfAway> = Context::new(160);

    // 1. AİLE A (Euler-Pfaff İkizleri):
    // a_n = 6n, b_n = -(2n - 1)(4n - (2k + 1))  (k = 0..6)
    let mut family_a_results = Vec::new();
    for k in 0..=6i64 {
        let a_fn = |n: i64| 6 * n;
        let b_fn = move |n: i64| -(2 * n - 1) * (4 * n - (2 * k + 1));

        let a_str = "6n".to_string();
        let b_str = format!("-(2n - 1) * (4n - {})", 2 * k + 1);

        if let Some(v_150) = evaluate_custom_pcf(a_fn, b_fn, 150, 160) {
            let v_f64 = f64::from_str(&v_150.to_string()[..20]).unwrap_or(0.0);
            let matches = find_matching_targets(v_f64, targets, 1.0e-5);

            let mut best_match = None;
            let mut best_diff = DBig::from(999999);

            for target in matches {
                let base_str = match target.base {
                    BaseConstant::Lemniscate => LEMNISCATE_STR,
                    BaseConstant::GaussAgm => GAUSS_AGM_STR,
                    BaseConstant::EllipticK => ELLIPTIC_K_STR,
                    BaseConstant::QuadraticRatio => QUADRATIC_RATIO_STR,
                };
                if let Ok(base_exact) = DBig::from_str(base_str) {
                    let (m1, m2, m3, m4) = target.m;
                    let dm1 = DBig::from(m1);
                    let dm2 = DBig::from(m2);
                    let dm3 = DBig::from(m3);
                    let dm4 = DBig::from(m4);

                    let num = ctx.mul(&dm1.repr(), &base_exact.repr()).value();
                    let num = ctx.add(&num.repr(), &dm2.repr()).value();
                    let den = ctx.mul(&dm3.repr(), &base_exact.repr()).value();
                    let den = ctx.add(&den.repr(), &dm4.repr()).value();

                    if den != DBig::from(0) {
                        let t_exact = ctx.div(&num.repr(), &den.repr()).value();
                        let diff = (&v_150 - &t_exact).abs();
                        if diff < best_diff {
                            best_diff = diff;
                            best_match = Some(target);
                        }
                    }
                }
            }

            let val_str = v_150.to_string();
            let display_val = if val_str.len() > 24 {
                format!("{}...", &val_str[..24])
            } else {
                val_str
            };

            let err_str = format_scientific(&best_diff);
            let matching_digits = parse_err_exp(&err_str) as usize;

            let (target_name, target_formula) = if let Some(m) = best_match {
                (m.name.clone(), m.formula.clone())
            } else {
                ("Bilinmiyor".to_string(), "-".to_string())
            };

            family_a_results.push(ParametricResult {
                k,
                a_str,
                b_str,
                value_str: display_val,
                target_name,
                target_formula,
                error_str: err_str,
                matching_digits,
            });
        }
    }

    // 2. AİLE B (Gauss-Pfaff Kaydırmaları):
    // a_n = 4n + (2k - 3), b_n = -n(2n - (2k - 1))  (k = 0..6)
    let mut family_b_results = Vec::new();
    for k in 0..=6i64 {
        let a_fn = move |n: i64| 4 * n + (2 * k - 3);
        let b_fn = move |n: i64| -n * (2 * n - (2 * k - 1));

        let a_str = if 2 * k - 3 > 0 {
            format!("4n + {}", 2 * k - 3)
        } else if 2 * k - 3 < 0 {
            format!("4n - {}", 3 - 2 * k)
        } else {
            "4n".to_string()
        };

        let b_str = if 2 * k - 1 > 0 {
            format!("-n * (2n - {})", 2 * k - 1)
        } else if 2 * k - 1 < 0 {
            format!("-n * (2n + {})", 1 - 2 * k)
        } else {
            "-n * 2n".to_string()
        };

        if let Some(v_150) = evaluate_custom_pcf(a_fn, b_fn, 150, 160) {
            let v_f64 = f64::from_str(&v_150.to_string()[..20]).unwrap_or(0.0);
            let matches = find_matching_targets(v_f64, targets, 1.0e-5);

            let mut best_match = None;
            let mut best_diff = DBig::from(999999);

            for target in matches {
                let base_str = match target.base {
                    BaseConstant::Lemniscate => LEMNISCATE_STR,
                    BaseConstant::GaussAgm => GAUSS_AGM_STR,
                    BaseConstant::EllipticK => ELLIPTIC_K_STR,
                    BaseConstant::QuadraticRatio => QUADRATIC_RATIO_STR,
                };
                if let Ok(base_exact) = DBig::from_str(base_str) {
                    let (m1, m2, m3, m4) = target.m;
                    let dm1 = DBig::from(m1);
                    let dm2 = DBig::from(m2);
                    let dm3 = DBig::from(m3);
                    let dm4 = DBig::from(m4);

                    let num = ctx.mul(&dm1.repr(), &base_exact.repr()).value();
                    let num = ctx.add(&num.repr(), &dm2.repr()).value();
                    let den = ctx.mul(&dm3.repr(), &base_exact.repr()).value();
                    let den = ctx.add(&den.repr(), &dm4.repr()).value();

                    if den != DBig::from(0) {
                        let t_exact = ctx.div(&num.repr(), &den.repr()).value();
                        let diff = (&v_150 - &t_exact).abs();
                        if diff < best_diff {
                            best_diff = diff;
                            best_match = Some(target);
                        }
                    }
                }
            }

            let val_str = v_150.to_string();
            let display_val = if val_str.len() > 24 {
                format!("{}...", &val_str[..24])
            } else {
                val_str
            };

            let err_str = format_scientific(&best_diff);
            let matching_digits = parse_err_exp(&err_str) as usize;

            let (target_name, target_formula) = if let Some(m) = best_match {
                (m.name.clone(), m.formula.clone())
            } else {
                ("Bilinmiyor".to_string(), "-".to_string())
            };

            family_b_results.push(ParametricResult {
                k,
                a_str,
                b_str,
                value_str: display_val,
                target_name,
                target_formula,
                error_str: err_str,
                matching_digits,
            });
        }
    }

    // 3. AİLE B2 (Discovered Gauss-Pfaff Alternating Parity Shift):
    // a_n = 4n + (2m - 3), b_n = -n * (2n - (5 - 2m))  (m = 0..4)
    let mut family_b2_results = Vec::new();
    for m in 0..=4i64 {
        let a_fn = move |n: i64| 4 * n + (2 * m - 3);
        let b_fn = move |n: i64| -n * (2 * n - (5 - 2 * m));

        let a_str = if 2 * m - 3 > 0 {
            format!("4n + {}", 2 * m - 3)
        } else if 2 * m - 3 < 0 {
            format!("4n - {}", 3 - 2 * m)
        } else {
            "4n".to_string()
        };

        let b_str = if 5 - 2 * m > 0 {
            format!("-n * (2n - {})", 5 - 2 * m)
        } else if 5 - 2 * m < 0 {
            format!("-n * (2n + {})", 2 * m - 5)
        } else {
            "-n * 2n".to_string()
        };

        if let Some(v_150) = evaluate_custom_pcf(a_fn, b_fn, 150, 160) {
            let v_f64 = f64::from_str(&v_150.to_string()[..20]).unwrap_or(0.0);
            let matches = find_matching_targets(v_f64, targets, 1.0e-5);

            let mut best_match = None;
            let mut best_diff = DBig::from(999999);

            for target in matches {
                let base_str = match target.base {
                    BaseConstant::Lemniscate => LEMNISCATE_STR,
                    BaseConstant::GaussAgm => GAUSS_AGM_STR,
                    BaseConstant::EllipticK => ELLIPTIC_K_STR,
                    BaseConstant::QuadraticRatio => QUADRATIC_RATIO_STR,
                };
                if let Ok(base_exact) = DBig::from_str(base_str) {
                    let (m1, m2, m3, m4) = target.m;
                    let dm1 = DBig::from(m1);
                    let dm2 = DBig::from(m2);
                    let dm3 = DBig::from(m3);
                    let dm4 = DBig::from(m4);

                    let num = ctx.mul(&dm1.repr(), &base_exact.repr()).value();
                    let num = ctx.add(&num.repr(), &dm2.repr()).value();
                    let den = ctx.mul(&dm3.repr(), &base_exact.repr()).value();
                    let den = ctx.add(&den.repr(), &dm4.repr()).value();

                    if den != DBig::from(0) {
                        let t_exact = ctx.div(&num.repr(), &den.repr()).value();
                        let diff = (&v_150 - &t_exact).abs();
                        if diff < best_diff {
                            best_diff = diff;
                            best_match = Some(target);
                        }
                    }
                }
            }

            let val_str = v_150.to_string();
            let display_val = if val_str.len() > 24 {
                format!("{}...", &val_str[..24])
            } else {
                val_str
            };

            let err_str = format_scientific(&best_diff);
            let matching_digits = parse_err_exp(&err_str) as usize;

            let (target_name, target_formula) = if let Some(target) = best_match {
                (target.name.clone(), target.formula.clone())
            } else {
                ("Bilinmiyor".to_string(), "-".to_string())
            };

            family_b2_results.push(ParametricResult {
                k: m,
                a_str,
                b_str,
                value_str: display_val,
                target_name,
                target_formula,
                error_str: err_str,
                matching_digits,
            });
        }
    }

    (family_a_results, family_b_results, family_b2_results)
}
