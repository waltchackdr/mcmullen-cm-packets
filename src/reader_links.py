"""Local navigation checks for maintained Markdown, without fetching remote URLs."""
import collections
from html import unescape
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


def markdown_anchors(text):
    anchors = set(re.findall(r'<(?:a|h[1-6])\b[^>]*\bid=["\']([^"\']+)', text))
    counts = collections.Counter()
    fenced = False
    for line in text.splitlines():
        if line.startswith(("```", "~~~")):
            fenced = not fenced
        if fenced:
            continue
        match = re.match(r"^#{1,6}\s+(.+?)(?:\s+#+)?$", line)
        if not match:
            continue
        title = re.sub(r"<[^>]*>", "", unescape(match[1])).lower()
        slug = re.sub(r"[^\w\- ]", "", title).replace(" ", "-")
        number = counts[slug]
        counts[slug] += 1
        anchors.add(slug + (f"-{number}" if number else ""))
    return anchors


def check_links(path, root, anchor_cache=None):
    """Check inline Markdown destinations and Markdown fragments within root.

    Remote URLs are not fetched. This intentionally does not claim to parse all
    Markdown syntax or validate fragments in PDFs, HTML or source-code viewers.
    """
    path, root = Path(path).resolve(), Path(root).resolve()
    cache = anchor_cache if anchor_cache is not None else {}
    errors = []
    links = fragments = 0
    for value in re.findall(r"(?<!!)\[[^\]]+\]\(([^)]+)\)", path.read_text()):
        target = value.strip("<>")
        parsed = urlsplit(target)
        if parsed.scheme or parsed.netloc:
            continue
        linked = (path.parent / unquote(parsed.path)).resolve() if parsed.path else path
        links += 1
        if not linked.is_relative_to(root) or not linked.exists():
            errors.append(f"Broken local link in {path.relative_to(root)}: {value}")
        elif parsed.fragment and linked.suffix.lower() == ".md":
            fragments += 1
            if linked not in cache:
                cache[linked] = markdown_anchors(linked.read_text())
            if unquote(parsed.fragment) not in cache[linked]:
                errors.append(f"Missing Markdown anchor in {path.relative_to(root)}: {value}")
    return links, fragments, errors
