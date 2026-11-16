"""Assemble the registered Julia package (Cantera.jl) from interfaces/julia.

The package differs from the source tree in two ways: it has no `src/generated`,
since the bindings, the library and the data files come from `Cantera_jll`, and
its Project.toml carries the Cantera version and depends on exactly that version
of `Cantera_jll`. The license is copied in, so that it has a single source.
"""

from pathlib import Path
import argparse
import shutil
import sys
import tomllib

# UUID that BinaryBuilder assigns to Cantera_jll (derived from its name).
CANTERA_JLL_UUID = "422db232-088f-5434-8d8c-b12d8b2bddc0"


def _toml_value(value) -> str:
    if isinstance(value, str):
        return f'"{value}"'
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    raise TypeError(f"unsupported TOML value {value!r}")


def project_toml(source: Path, version: str, jll_version: str) -> str:
    project = tomllib.loads(source.read_text())
    project["version"] = version
    project.setdefault("deps", {})["Cantera_jll"] = CANTERA_JLL_UUID
    project.setdefault("compat", {})["Cantera_jll"] = f"={jll_version}"

    lines = [f"{key} = {_toml_value(project[key])}"
             for key in ("name", "uuid", "authors", "version")]
    for table in ("deps", "compat", "extras", "targets"):
        if table in project:
            lines += ["", f"[{table}]"]
            lines += [f"{key} = {_toml_value(value)}"
                      for key, value in sorted(project[table].items())]
    return "\n".join(lines) + "\n"


def main(source_directory: Path, target_directory: Path, version: str,
         jll_version: str) -> Path:
    julia_source = source_directory / "interfaces" / "julia"
    if target_directory.exists():
        shutil.rmtree(target_directory)
    shutil.copytree(julia_source / "src", target_directory / "src",
                    ignore=shutil.ignore_patterns("generated"))
    shutil.copytree(julia_source / "test", target_directory / "test")
    shutil.copy2(julia_source / "README.md", target_directory / "README.md")
    shutil.copy2(source_directory / "License.txt", target_directory / "LICENSE")
    (target_directory / "Project.toml").write_text(
        project_toml(julia_source / "Project.toml", version, jll_version))
    return target_directory


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_directory")
    parser.add_argument("target_directory")
    parser.add_argument("version", help="version of Cantera.jl, e.g. 4.0.0")
    parser.add_argument("--jll-version", help="Cantera_jll version (default: same)")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    print(main(Path(args.source_directory), Path(args.target_directory),
               args.version, args.jll_version or args.version))
