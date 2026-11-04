"""Bounded reader tooling. This module never executes mathematical certificates."""
from pathlib import Path
import hashlib
import re

ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / "papers"
DOCUMENTS = {
    "packet": "packet-and-principal-models",
    "quotients": "quotients-and-descent",
    "hodge": "hodge-classes",
    "supplement": "arithmetic-and-computation",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tex_inputs(path, base=PAPERS, stack=()):
    """Expand the actual input graph, rejecting cycles and escaping paths."""
    path = path.resolve()
    if not path.is_relative_to(base.resolve()) or path in stack:
        raise ValueError(f"Invalid TeX input graph at {path}")
    source = path.read_text()

    def expand(match):
        target = base / match[1]
        if not target.suffix:
            target = target.with_suffix(".tex")
        return tex_inputs(target, base, (*stack, path))

    return re.sub(r"\\input\{([^}]+)\}", expand, source)


def auxiliary_numbers(path):
    """Read the first balanced group of each LaTeX newlabel value."""
    source = Path(path).read_text()
    result = {}
    for match in re.finditer(r"\\newlabel\{([^}]+)\}\{\{", source):
        label = match[1]
        start = end = match.end()
        depth = 1
        while depth and end < len(source):
            if source[end] == "{":
                depth += 1
            elif source[end] == "}":
                depth -= 1
            end += 1
        if depth:
            raise ValueError(f"Unbalanced auxiliary label: {label}")
        if "@cref" not in label:
            result[label] = source[start:end - 1]
    return result
