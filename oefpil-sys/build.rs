//! Build script for C/Fortran library OEFPIL.

use cc::Build;
use std::path::Path;

#[allow(clippy::too_many_lines)]
fn main() {
    let src = Path::new("src");
    let changed = src.to_str().unwrap();
    println!("cargo::rerun-if-changed={changed}");
    if cfg!(feature = "built-in") {
        fn gcc_library_path() -> Option<String> {
            Build::new()
                .compiler("gcc")
                .get_compiler()
                .to_command()
                .arg("-print-libgcc-file-name")
                .output()
                .ok()
                .and_then(|output| String::from_utf8(output.stdout).ok())
                .and_then(|utf8| {
                    Path::new(utf8.trim_end())
                        .parent()
                        .and_then(Path::to_str)
                        .map(str::to_string)
                })
        }
        let lapack_path = src.join("lapack");
        let blas_path = lapack_path.join("blas");
        #[rustfmt::skip]
		let blas_file = [
			"dcopy.f",
			"ddot.f",
			"dgemm.f",
			"dgemv.f",
			"dger.f",
			"dnrm2.f90",
			"dscal.f",
			"dsymm.f",
			"dsymv.f",
			"dsyrk.f",
			"dtrmm.f",
			"dtrmv.f",
			"dtrsm.f",
			"dtrsv.f",
			"lsame.f",
		]
		.map(|f| blas_path.join(f));
        #[rustfmt::skip]
		let lapack_file = [
			"dgemqr.f",
			"dgemqrt.f",
			"dgeqr.f",
			"dgeqrt.f",
			"dgeqrt2.f",
			"dgeqrt3.f",
			"disnan.f",
			"dlaisnan.f",
			"dlamch.f",
			"dlamtsqr.f",
			"dlapy2.f",
			"dlarfb.f",
			"dlarfg.f",
			"dlatsqr.f",
			"dlauu2.f",
			"dlauum.f",
			"dpotrf.f",
			"dpotrf2.f",
			"dpotri.f",
			"dtpmqrt.f",
			"dtpqrt.f",
			"dtpqrt2.f",
			"dtprfb.f",
			"dtrti2.f",
			"dtrtri.f",
			"ieeeck.f",
			"ilaenv.f",
			"iparmq.f",
			"xerbla.f",
		]
		.map(|f| lapack_path.join(f));
        Build::new()
            .compiler("gcc")
            .files(blas_file)
            .flag("-Wno-compare-reals")
            .flag("-Wno-function-elimination")
            .flag_if_supported("-Wno-maybe-uninitialized")
            .warnings_into_errors(true)
            .compile("oefpil-blas");
        Build::new()
            .compiler("gcc")
            .files(lapack_file)
            .flag("-Wno-compare-reals")
            .flag("-Wno-unused-dummy-argument")
            .flag("-Wno-function-elimination")
            .flag_if_supported("-Wno-maybe-uninitialized")
            .warnings_into_errors(true)
            .compile("oefpil-lapack");
        println!("cargo::rustc-link-lib=static=oefpil-blas");
        println!("cargo::rustc-link-lib=static=oefpil-lapack");
        if let Some(path) = gcc_library_path() {
            println!("cargo::rustc-link-search={path}");
            println!("cargo::rustc-link-lib=static=gfortran");
            println!("cargo::rustc-link-lib=static=quadmath");
        } else {
            println!("cargo::rustc-link-lib=gfortran");
        }
    } else {
        println!("cargo::rustc-link-lib=lapack");
        println!("cargo::rustc-link-lib=blas");
    }
    let oefpil_path = src.join("oefpil");
    #[rustfmt::skip]
	let oefpil_file = [
		"blaslapack.c",
		"oefpil.c",
		"oefpil_tcm.c",
		"tiled_la.c",
	]
	.map(|c| oefpil_path.join(c));
    let lib_path = src;
    #[rustfmt::skip]
	let lib_file = [
		"lib.c",
	]
	.map(|c| lib_path.join(c));
    Build::new()
        .include(&oefpil_path)
        .files(oefpil_file)
        .include(lib_path)
        .files(&lib_file)
        .flag("-Wno-unused-variable")
        .flag("-Wno-unused-parameter")
        .flag_if_supported("-Wno-maybe-uninitialized")
        .warnings_into_errors(true)
        .compile("oefpil");
    #[cfg(feature = "bindgen")]
    {
        use bindgen::{Builder, callbacks::ParseCallbacks};
        use std::{env::var, path::PathBuf};
        #[derive(Debug)]
        struct RenameItem<'a> {
            old: &'a str,
            new: &'a str,
        }
        impl ParseCallbacks for RenameItem<'_> {
            fn item_name(&self, old: &str) -> Option<String> {
                (self.old == old).then(|| self.new.into())
            }
        }
        let out_path = PathBuf::from(var("OUT_DIR").unwrap());
        Builder::default()
            .header(lib_path.join("lib.h").to_str().unwrap())
            .allowlist_type("oefpil.*")
            .allowlist_function("oefpil.*")
            .allowlist_function("std.*")
            .blocklist_type("FILE")
            .opaque_type("FILE")
            .parse_callbacks(Box::new(RenameItem {
                old: "fitfunc_oefpil",
                new: "Evaluate",
            }))
            .generate()
            .expect("Cannot generate bindings")
            .write_to_file(out_path.join("bindings.rs"))
            .expect("Cannot write bindings");
    }
}
