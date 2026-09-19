"""Load and validate the master parameter file.

Every parameter in config.yaml is a mapping with at least ``value`` and
``status``.  The loader exists mainly to make two things impossible to do
by accident:

  1. use a PLACEHOLDER parameter whose value is still ``null``, and
  2. lose track of which numbers are measured, cited or invented.

``audit()`` is the honesty report: it lists what still has to be measured.
"""

from __future__ import annotations

import os
from typing import Any

import yaml

VALID_STATUS = {
    "MEASURED",
    "PLACEHOLDER",
    "LITERATURE",
    "ASSUMPTION",
    "DERIVED",
    "NUMERICAL",
}

DEFAULT_CONFIG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "config",
    "config.yaml",
)


def load(path: str = DEFAULT_CONFIG) -> dict:
    """Parse the YAML config."""
    with open(path) as fh:
        return yaml.safe_load(fh)


def _is_param(node: Any) -> bool:
    return isinstance(node, dict) and "status" in node


def walk(cfg: dict):
    """Yield ``(dotted_name, param_dict)`` for every parameter entry."""
    for section, body in cfg.items():
        if not isinstance(body, dict):
            continue
        for name, node in body.items():
            if _is_param(node):
                yield f"{section}.{name}", node


def get(cfg: dict, dotted: str) -> Any:
    """Fetch a parameter value, refusing to hand back an unset placeholder."""
    section, _, name = dotted.partition(".")
    try:
        node = cfg[section][name]
    except KeyError as exc:
        raise KeyError(f"no such parameter: {dotted}") from exc
    if not _is_param(node):
        raise KeyError(f"{dotted} is not a parameter entry")
    if node.get("value") is None:
        raise ValueError(
            f"{dotted} is an unset {node['status']}. "
            "Supply a measured value before using it."
        )
    return node["value"]


def validate(cfg: dict) -> list[str]:
    """Return a list of structural problems. Empty list means the file is sane."""
    problems: list[str] = []
    for dotted, node in walk(cfg):
        status = node["status"]
        if status not in VALID_STATUS:
            problems.append(f"{dotted}: unknown status {status!r}")
        if status in {"PLACEHOLDER", "ASSUMPTION", "LITERATURE"} and not node.get(
            "source"
        ):
            problems.append(f"{dotted}: status {status} requires a 'source'")
        value, rng = node.get("value"), node.get("range")
        if rng is not None:
            if len(rng) != 2 or rng[0] > rng[1]:
                problems.append(f"{dotted}: malformed range {rng!r}")
            elif isinstance(value, (int, float)) and not rng[0] <= value <= rng[1]:
                problems.append(f"{dotted}: default {value} outside range {rng}")
    return problems


def audit(cfg: dict) -> dict[str, list[str]]:
    """Group parameter names by status, for the 'what do I still need' report."""
    out: dict[str, list[str]] = {}
    for dotted, node in walk(cfg):
        out.setdefault(node["status"], []).append(dotted)
    return out


def unset_placeholders(cfg: dict) -> list[str]:
    """PLACEHOLDER parameters that still have no value at all."""
    return [
        dotted
        for dotted, node in walk(cfg)
        if node["status"] == "PLACEHOLDER" and node.get("value") is None
    ]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    cfg = load()
    problems = validate(cfg)
    print(f"config: {sum(1 for _ in walk(cfg))} parameters, {len(problems)} problems")
    for p in problems:
        print("  PROBLEM:", p)
    print()
    for status, names in sorted(audit(cfg).items()):
        print(f"{status} ({len(names)})")
        for n in names:
            print("   ", n)
    print()
    print("STILL UNMEASURED (value is null):")
    for n in unset_placeholders(cfg):
        print("   ", n)
