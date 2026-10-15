use catalan_pcf_miner::pcf::{evaluate_wallis_35, prefilter, PcfCandidate, PcfFamily};
use catalan_pcf_miner::verify::{
    build_target_table, find_matching_targets, save_identity, verify_stage3, TargetEntry,
};
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashSet;
use std::env;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct SearchConfig {
    pub threads: usize,
    pub mobius_depth: i64,
    pub family: String,
    pub tolerance: f64,
    pub output_file: String,

    // Aile 1 sınırları
    pub f1_c1_max: i64,
    pub f1_c2_min: i64,
    pub f1_c2_max: i64,
    pub f1_c3_max: i64,
    pub f1_q_min: i64,
    pub f1_q_max: i64,

    // Aile 2 sınırları
    pub f2_c1_max: i64,
    pub f2_c2_min: i64,
    pub f2_c2_max: i64,
    pub f2_c3_min: i64,
    pub f2_c3_max: i64,
    pub f2_c4_max: i64,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            threads: 7,
            mobius_depth: 6,           // [-6, 6] Möbius derinliği
            family: "all".to_string(), // 1, 2, all
            tolerance: 1.0e-7,
            output_file: "discovered_lemniscate_identities.txt".to_string(),

            // Aile 1: Lineer Gauss-Pfaff
            f1_c1_max: 16,
            f1_c2_min: -12,
            f1_c2_max: 12,
            f1_c3_max: 8,
            f1_q_min: -5,
            f1_q_max: 5,

            // Aile 2: Kuadratik Contiguous
            f2_c1_max: 16,
            f2_c2_min: -16,
            f2_c2_max: 16,
            f2_c3_min: -16,
            f2_c3_max: 16,
            f2_c4_max: 6,
        }
    }
}

fn print_help() {
    println!(
        r#"
================================================================================
 Mod-4 Hipergeometrik Lemniscate ve Eliptik Sabitler PCF Arama Motoru
================================================================================
Kullanım:
  catalan-pcf-miner [SEÇENEKLER]

Hedef Sabitler:
  1. Gauss Lemniscate Sabiti (ϖ ≈ 2.62205755...)
  2. Gauss AGM Sabiti (pi / ϖ ≈ 1.19814023...)
  3. İkinci Lemniscate / Eliptik Sabit K(1/√2) (ϖ / √2 ≈ 1.85407467...)
  4. Karesel Lemniscate Oranı (ϖ² / pi ≈ 2.18843961...)

Aileler:
  --family <all|1|2>
    all: Her iki aile (1 ve 2)
    1:   Aile 1: Asimetrik Çarpanlı Lineer Gauss-Pfaff Ailesi
    2:   Aile 2: Kuadratik a_n ve 4. Derece Mod-4 Çarpanlı b_n (Contiguous)

Seçenekler:
  -t, --threads <N>         Rayon thread sayısı (varsayılan: 7, max: 8)
  -m, --mobius-depth <N>    Möbius derinliği (varsayılan: 6)
  -tol, --tolerance <F>     Aşama 2 toleransı (varsayılan: 1e-7)
  -o, --output <DOSYA>      Kayıt dosyası (varsayılan: discovered_lemniscate_identities.txt)
  -h, --help                Bu yardım menüsünü görüntüler
"#
    );
}

fn parse_args() -> SearchConfig {
    let mut config = SearchConfig::default();
    let args: Vec<String> = env::args().collect();
    let mut i = 1;

    while i < args.len() {
        match args[i].as_str() {
            "-t" | "--threads" => {
                if i + 1 < args.len() {
                    config.threads = args[i + 1].parse().unwrap_or(config.threads);
                    i += 1;
                }
            }
            "-m" | "--mobius-depth" => {
                if i + 1 < args.len() {
                    config.mobius_depth = args[i + 1].parse().unwrap_or(config.mobius_depth);
                    i += 1;
                }
            }
            "--family" => {
                if i + 1 < args.len() {
                    config.family = args[i + 1].to_lowercase();
                    i += 1;
                }
            }
            "-tol" | "--tolerance" => {
                if i + 1 < args.len() {
                    config.tolerance = args[i + 1].parse().unwrap_or(config.tolerance);
                    i += 1;
                }
            }
            "-o" | "--output" => {
                if i + 1 < args.len() {
                    config.output_file = args[i + 1].clone();
                    i += 1;
                }
            }
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            _ => {}
        }
        i += 1;
    }

    config
}

#[inline(always)]
fn evaluate_and_verify(
    candidate: &PcfCandidate,
    targets: &[TargetEntry],
    tolerance: f64,
    seen_identities: &Mutex<HashSet<(String, String)>>,
    output_file: &str,
    pb: &ProgressBar,
    float_passed: &AtomicUsize,
    identities_found: &AtomicUsize,
) {
    if let Some(val) = evaluate_wallis_35(candidate) {
        let matching = find_matching_targets(val, targets, tolerance);
        if !matching.is_empty() {
            float_passed.fetch_add(1, Ordering::Relaxed);
            for target in matching {
                if let Some(res) = verify_stage3(candidate, target) {
                    let mut seen = seen_identities.lock().unwrap();
                    let key = (res.a_repr.clone(), res.b_repr.clone());
                    if seen.insert(key) {
                        identities_found.fetch_add(1, Ordering::Relaxed);

                        pb.suspend(|| {
                            println!("\n[BULDUM!]");
                            println!("Hedef     : {}", res.target_name);
                            println!("a_n       : {}", res.a_repr);
                            println!("b_n       : {}", res.b_repr);
                            println!("Değer     : {}", res.value_str);
                            println!("Hata Payı : {}", res.error_str);
                            println!("Adım      : {} iterasyon\n", res.steps);
                        });

                        let _ = save_identity(&res, output_file);
                    }
                }
            }
        }
    }
}

fn main() {
    let config = parse_args();

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .expect("Rayon thread pool başlatılamadı!");

    println!("================================================================================");
    println!(" Mod-4 Hipergeometrik Lemniscate ve Eliptik Sabitler PCF Arama Motoru");
    println!(" Hedef Mimari: Intel Core i3-N305 (8 Gracemont E-Core, AVX2, 6MB L3)");
    println!("================================================================================");
    println!("Aktif İş Parçacığı (Rayon Threads) : {}", config.threads);
    println!("Möbius Dönüşüm Derinliği            : [-{}, {}]", config.mobius_depth, config.mobius_depth);
    println!("Taranacak Aile(ler)                 : {}", config.family);
    println!("Aşama 2 Toleransı                   : {:e}", config.tolerance);
    println!("Kayıt Dosyası                       : {}", config.output_file);

    let _ = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .open(&config.output_file);

    // Hedef Tablosunu Oluştur (4 temel sabit ve [-6, 6] Möbius dönüşümleri)
    let targets = build_target_table(config.mobius_depth);
    println!("Hedef Möbius Tablosu Boyutu         : {} sabit", targets.len());

    let run_f1 = config.family == "all" || config.family == "1";
    let run_f2 = config.family == "all" || config.family == "2";

    // Sabit parametre kümeleri
    // Aile 1 için (p1, p2) çiftleri
    let p_pairs_f1: [(i64, i64); 6] = [(1, 1), (1, 2), (1, 4), (2, 2), (2, 4), (4, 4)];

    // Aile 2 için (p1, p2) çiftleri ve q değerleri (tek sayılar: [-3, -1, 1, 3])
    let p_pairs_f2: [(i64, i64); 4] = [(1, 2), (2, 2), (2, 4), (4, 4)];
    let q_vals_f2: [i64; 4] = [-3, -1, 1, 3];

    // Kombinasyon hesaplama
    let mut f1_combinations = 0u64;
    if run_f1 {
        let c1_count = config.f1_c1_max as u64;
        let c2_count = (config.f1_c2_max - config.f1_c2_min + 1) as u64;
        let c3_count = config.f1_c3_max as u64;
        let q_len = (config.f1_q_max - config.f1_q_min + 1) as u64;

        let mut q_pairs_count = 0u64;
        for &(p1, p2) in &p_pairs_f1 {
            if p1 == p2 {
                q_pairs_count += q_len * (q_len + 1) / 2;
            } else {
                q_pairs_count += q_len * q_len;
            }
        }
        f1_combinations = c1_count * c2_count * c3_count * q_pairs_count;
    }

    let mut f2_combinations = 0u64;
    if run_f2 {
        let c1_count = config.f2_c1_max as u64;
        let c2_count = (config.f2_c2_max - config.f2_c2_min + 1) as u64;
        let c3_count = (config.f2_c3_max - config.f2_c3_min + 1) as u64;
        let c4_count = config.f2_c4_max as u64;
        let k_count = 3u64; // k in [0, 2]

        let mut q_pairs_f2_count = 0u64;
        for &(p1, p2) in &p_pairs_f2 {
            if p1 == p2 {
                q_pairs_f2_count += 4 * 5 / 2; // 10
            } else {
                q_pairs_f2_count += 4 * 4; // 16
            }
        }
        f2_combinations = c1_count * c2_count * c3_count * c4_count * k_count * q_pairs_f2_count;
    }

    let total_combinations = f1_combinations + f2_combinations;
    if run_f1 {
        println!("Aile 1 (Lineer Gauss-Pfaff)         : {} kombinasyon", f1_combinations);
    }
    if run_f2 {
        println!("Aile 2 (Kuadratik Contiguous)       : {} kombinasyon", f2_combinations);
    }
    println!("Toplam Taranacak Kombinasyon        : {}", total_combinations);
    println!("================================================================================\n");

    let pb = ProgressBar::new(total_combinations);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) | Hız: {per_sec} | Kalan: {eta}")
            .unwrap()
            .progress_chars("#>-"),
    );

    let prefilter_passed = AtomicUsize::new(0);
    let float_passed = AtomicUsize::new(0);
    let identities_found = AtomicUsize::new(0);
    let seen_identities = Mutex::new(HashSet::new());

    let start_time = Instant::now();

    use rayon::prelude::*;

    // 1. AİLE 1: ASİMETRİK ÇARPANLI LİNEER GAUSS-PFAFF AİLESİ
    // a_n = c1*n + c2
    // b_n = -c3 * (p1*n + q1) * (p2*n + q2)
    if run_f1 {
        let mut tasks_f1 = Vec::new();
        for c1 in 1..=config.f1_c1_max {
            for c2 in config.f1_c2_min..=config.f1_c2_max {
                if c1 + c2 <= 0 {
                    continue; // a_1 pozitif olmalı
                }
                tasks_f1.push((c1, c2));
            }
        }

        tasks_f1.par_iter().for_each(|&(c1, c2)| {
            let mut local_pre_pass = 0usize;
            let mut local_scanned = 0u64;
            let c1_sq = c1 * c1;

            for &(p1, p2) in &p_pairs_f1 {
                for c3 in 1..=config.f1_c3_max {
                    // Darboux limit testi: 4 * c3 * p1 * p2 <= c1^2
                    let factor_4 = 4 * c3 * p1 * p2;
                    let darboux_ok = factor_4 <= c1_sq;

                    for q1 in config.f1_q_min..=config.f1_q_max {
                        let q2_start = if p1 == p2 { q1 } else { config.f1_q_min };
                        for q2 in q2_start..=config.f1_q_max {
                            local_scanned += 1;

                            if !darboux_ok {
                                continue;
                            }

                            let cand = PcfCandidate {
                                family: PcfFamily::LinearGaussPfaff,
                                c1,
                                c2,
                                c3: 0,
                                c_factor: c3,
                                k: 0,
                                p1,
                                q1,
                                p2,
                                q2,
                            };

                            if prefilter(&cand) {
                                local_pre_pass += 1;
                                evaluate_and_verify(
                                    &cand,
                                    &targets,
                                    config.tolerance,
                                    &seen_identities,
                                    &config.output_file,
                                    &pb,
                                    &float_passed,
                                    &identities_found,
                                );
                            }
                        }
                    }
                }
            }

            prefilter_passed.fetch_add(local_pre_pass, Ordering::Relaxed);
            pb.inc(local_scanned);
        });
    }

    // 2. AİLE 2: KUADRATİK a_n VE 4. DERECE MOD-4 ÇARPANLI b_n (CONTIGUOUS)
    // a_n = c1*n^2 + c2*n + c3
    // b_n = -c4 * n * (n + k) * (p1*n + q1) * (p2*n + q2)
    if run_f2 {
        let mut tasks_f2 = Vec::new();
        for c1 in 1..=config.f2_c1_max {
            for c2 in config.f2_c2_min..=config.f2_c2_max {
                tasks_f2.push((c1, c2));
            }
        }

        tasks_f2.par_iter().for_each(|&(c1, c2)| {
            let mut local_pre_pass = 0usize;
            let mut local_scanned = 0u64;
            let c1_sq = c1 * c1;

            for c3 in config.f2_c3_min..=config.f2_c3_max {
                if c1 + c2 + c3 <= 0 {
                    // a_1 <= 0 olan durumlar
                    // Tüm alt kombinasyonları taranmış say
                    let skipped = (config.f2_c4_max as u64) * 3 * 52;
                    local_scanned += skipped;
                    continue;
                }

                for &(p1, p2) in &p_pairs_f2 {
                    for c4 in 1..=config.f2_c4_max {
                        // Darboux limit testi: 4 * c4 * p1 * p2 <= c1^2
                        let factor_4 = 4 * c4 * p1 * p2;
                        let darboux_ok = factor_4 <= c1_sq;

                        for k in 0..=2i64 {
                            for &q1 in &q_vals_f2 {
                                for &q2 in &q_vals_f2 {
                                    if p1 == p2 && q1 > q2 {
                                        continue;
                                    }

                                    local_scanned += 1;

                                    if !darboux_ok {
                                        continue;
                                    }

                                    let cand = PcfCandidate {
                                        family: PcfFamily::QuadraticContiguous,
                                        c1,
                                        c2,
                                        c3,
                                        c_factor: c4,
                                        k,
                                        p1,
                                        q1,
                                        p2,
                                        q2,
                                    };

                                    if prefilter(&cand) {
                                        local_pre_pass += 1;
                                        evaluate_and_verify(
                                            &cand,
                                            &targets,
                                            config.tolerance,
                                            &seen_identities,
                                            &config.output_file,
                                            &pb,
                                            &float_passed,
                                            &identities_found,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }

            prefilter_passed.fetch_add(local_pre_pass, Ordering::Relaxed);
            pb.inc(local_scanned);
        });
    }

    pb.finish();

    let duration = start_time.elapsed();
    let total_scanned = total_combinations;
    let pre_pass = prefilter_passed.load(Ordering::Relaxed);
    let flt_pass = float_passed.load(Ordering::Relaxed);
    let ids_found = identities_found.load(Ordering::Relaxed);

    let secs = duration.as_secs_f64();
    let throughput = if secs > 0.0 {
        total_scanned as f64 / secs
    } else {
        0.0
    };

    println!("\n================================================================================");
    println!(" TARAMA ÖZETİ");
    println!("================================================================================");
    println!("Toplam Taranan Kombinasyon : {}", total_scanned);
    println!(
        "Aşama 1 (Ön-Eleme) Geçen   : {} ({:.2}%)",
        pre_pass,
        if total_scanned > 0 {
            (pre_pass as f64 / total_scanned as f64) * 100.0
        } else {
            0.0
        }
    );
    println!("Aşama 2 (f64 Wallis) Geçen : {}", flt_pass);
    println!("Aşama 3 Doğrulanan Kimlik  : {}", ids_found);
    println!("Toplam Süre                : {:?}", duration);
    println!("Arama Hızı (Throughput)    : {:.0} kombinasyon/saniye", throughput);
    println!("Kayıt Dosyası              : {}", config.output_file);
    println!("================================================================================");
}
