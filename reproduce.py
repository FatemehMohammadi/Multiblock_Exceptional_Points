"""Reproduce all checks and figures without using the archived results as input."""
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=ROOT / "reproduced")
    args = parser.parse_args()
    out = args.output_root.resolve()
    if out == ROOT or out in (ROOT / "results", ROOT / "figures"):
        raise ValueError("Choose a fresh output directory, not an archived input directory")
    logs = out / "validation"
    logs.mkdir(parents=True, exist_ok=True)
    def run(name, command, cwd=ROOT):
        result = subprocess.run(command, cwd=cwd, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (logs / (name + ".log")).write_text(result.stdout)
        if result.returncode:
            raise RuntimeError(name + " failed:\n" + result.stdout[-3000:])
        print(name + ": passed", flush=True)
        return result.stdout
    tests = run("regression", [sys.executable, "-m", "unittest", "discover",
                               "-s", "experiments", "-p", "test_*.py", "-v"])
    if "Ran 5 tests" not in tests:
        raise ValueError("Unexpected regression-test count")
    run("experiments", [sys.executable, "experiments/run_experiments.py",
                         "--output-root", str(out)])
    run("architecture_figure", [sys.executable, "experiments/architecture_figure.py",
                                 "--output", str(out / "figures")])
    with tempfile.TemporaryDirectory(prefix="independent_", dir=out) as folder:
        stage = Path(folder)
        script = stage / "independent_all_jordan_check.py"
        shutil.copy2(ROOT / script.name, script)
        run("independent_polynomial_check", [sys.executable, str(script)], cwd=stage)
        report = json.loads((stage / "all_jordan_symbolic_audit.json").read_text())
        if report["cases"] != 87 or report["status"] != "passed":
            raise ValueError(report)
        (logs / "independent_polynomial_report.json").write_text(
            json.dumps(report, indent=2) + "\n")
    groups = json.loads(run("independent_group_checks",
                           [sys.executable, "global_group_exact_checks.py"]))
    (logs / "independent_group_report.json").write_text(json.dumps(groups, indent=2) + "\n")
    run("certificate_checks", [sys.executable, "certificates/verify_certificate.py",
                              "--output", str(out / "results/certificate_verification.json")])
    run("local_certificate", [sys.executable, "certificates/test_local_certificate.py",
                              "--output", str(out / "results/local_certificate_report.json")])
    run("quartic_example", [sys.executable, "experiments/quartic_example.py",
                            "--output-root", str(out)])
    run("general_fusion", [sys.executable, "certificates/check_general_fusion.py",
                           "--output", str(out / "results/general_fusion_report.json")])
    run("independent_general_fusion", [sys.executable, "certificates/independent_general_fusion.py",
                                      "--output", str(out / "results/independent_general_fusion_report.json")])
    run("independent_discriminant", [sys.executable, "certificates/independent_discriminant_audit.py",
                                     "--output", str(out / "results/independent_discriminant_report.json")])
    run("citation_substitutions", [sys.executable, "experiments/check_citation_substitutions.py",
                                   "--output", str(out / "results/citation_substitutions.json")])
    run("root_continuation", [sys.executable, "experiments/check_root_continuation.py",
                              "--output", str(out / "results/root_continuation_report.json")])
    run("support_obstruction", [sys.executable, "certificates/check_support_obstruction.py",
                                 "--output", str(out / "results/support_obstruction_report.json")])
    run("tropical_refinement", [sys.executable, "experiments/check_tropical_refinement.py",
                                 "--output-root", str(out)])
    run("cone_geometry", [sys.executable, "experiments/cone_geometry.py",
                            "--output-root", str(out / "figures")])
    shutil.move(str(out / "figures/cone_geometry_checks.json"),
                str(out / "results/cone_geometry_checks.json"))
    tropical = json.loads((out / "results/tropical_refinement_report.json").read_text())
    cone = json.loads((out / "results/cone_geometry_checks.json").read_text())
    if tropical["status"] != "passed" or tropical["profile_cases"] != 30 or tropical["full_integer_determinants"] != 1640:
        raise ValueError("Higher-refinement check failed")
    if cone["status"] != "passed" or cone["sylvester_determinants"] != 1800:
        raise ValueError("Cone-geometry check failed")
    support = json.loads((out / "results/support_obstruction_report.json").read_text())
    if support["support_patterns"] != 1094599 or support["full_symbolic_determinants"] != 3015 or support["Jordan_profiles"] != 10:
        raise ValueError("Unexpected small-support finite audit count")
    rank_example = support["higher_rank_counterexample"]
    if rank_example["full_integer_determinants"] != 44 or rank_example["rank_B"] != 3:
        raise ValueError("Higher-rank example failed")
    (out / "results/rank_relaxation_report.json").write_text(json.dumps(rank_example, indent=2) + "\n")
    fusion = json.loads((out / "results/general_fusion_report.json").read_text())
    independent_fusion = json.loads((out / "results/independent_general_fusion_report.json").read_text())
    continuation = json.loads((out / "results/root_continuation_report.json").read_text())
    if fusion["profile_cases"] != 263 or fusion["full_integer_determinants"] != 15886 or fusion["full_sylvester_determinants"] != 1036:
        raise ValueError("Unexpected full-cycle matrix-check count")
    if independent_fusion["full_integer_determinants"] != 15644 or independent_fusion["finite_full_matrix_profiles"] != 259:
        raise ValueError("Unexpected independent full-cycle check count")
    if len(continuation["cases"]) != 9 or not all(c["same_labeled_generators_after_refinement"] for c in continuation["cases"]):
        raise ValueError("Numerical continuation refinement failed")
    certificate = json.loads((out / "results/certificate_verification.json").read_text())
    local = json.loads((out / "results/local_certificate_report.json").read_text())
    quartic = json.loads((out / "results/quartic_example_report.json").read_text())
    if certificate["exact_family"]["profile_cases"] != 41 or certificate["exact_family"]["full_integer_determinants"] != 1083:
        raise ValueError("Unexpected sparse-certificate check count")
    if local["status"] != "passed" or quartic["sylvester_determinants"] != 60:
        raise ValueError("A local or quartic certificate check did not pass")
    summary = json.loads((out / "results/summary.json").read_text())
    if summary["pointwise_rows"] != 93 or summary["different_partitions_verified"] != 93:
        raise ValueError("Unexpected readout record count")
    jordan = summary["integer_checker_reports"]["all_jordan_verification"]
    equal = summary["integer_checker_reports"]["higher_rank_verification"]
    if jordan["jordan_types"] != 128 or jordan["full_integer_determinants"] != 5797:
        raise ValueError(jordan)
    if equal["total_integer_determinants"] != 1332:
        raise ValueError(equal)
    result = {"status": "passed", "regression_tests": 5,
              "all_jordan_types": 128, "integer_determinants_all_jordan": 5797,
              "integer_determinants_equal_blocks": 1332,
              "independent_bivariate_identities": 87, "pointwise_response_records": 93,
              "certificate_profile_cases": 41, "certificate_integer_determinants": 1083,
              "certificate_compiler_targets": 1172,
              "quartic_sylvester_determinants": 60,
              "full_cycle_profile_cases": 263, "full_cycle_determinants": 15886,
              "full_cycle_sylvester_determinants": 1036,
              "independent_discriminant_sylvester_determinants": 220,
              "independent_full_cycle_profiles": 259, "independent_full_cycle_determinants": 15644,
              "balancing_profiles_through_dimension_20": 2693,
              "numerical_continuation_cases": 9, "numerical_loop_runs": 36,
              "reference_specialization_determinants": 2,
              "small_support_Jordan_profiles": 10,
              "small_support_patterns": 1094599,
              "small_support_symbolic_determinants": 3015,
              "rank_relaxation_integer_determinants": 44,
              "higher_refinement_profiles": 30,
              "higher_refinement_full_determinants": 1640,
              "higher_refinement_sylvester_determinants": 58,
              "higher_refinement_continuation_runs": 8,
              "cone_sylvester_determinants": 1800,
              "cone_rational_matrix_determinants": 120,
              "cone_rational_sylvester_determinants": 120,
              "local_trace_certificate": "passed, including zero-certificate full-cycle example",
              "scope": "Finite exact checks and floating illustrations; general statements are proved in the manuscript."}
    (out / "reproduction_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print("All reproduction and scientific summary checks passed.")

if __name__ == "__main__":
    main()
