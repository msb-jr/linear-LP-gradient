"""Figure/table notes for outputs.

Each script calls `build_notes()` with a list of legends, producing
a plain text file with captions and footnotes for each output.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
import textwrap


def _ts_iso(timestamp: Optional[datetime] = None) -> str:
    ts = timestamp or datetime.now(timezone.utc)
    return ts.astimezone(timezone.utc).replace(microsecond=0).isoformat()


def _wrap(text: str, width: int = 78, indent: str = "") -> str:
    """Wrap text to width with optional indent for continuation lines."""
    lines = textwrap.wrap(text, width=width - len(indent))
    if not lines:
        return ""
    return lines[0] + "\n" + "\n".join(indent + line for line in lines[1:])


def build_notes(legends: List[Dict[str, Any]], title: Optional[str] = None) -> str:
    """Build a notes file from a list of figure/table legends.

    Parameters
    ----------
    legends : list of dict
        Each dict should have:
          - target: str  (e.g., "Figure 1", "Table 1", or filename)
          - caption: str (the legend/caption text)
          - footnotes: list[str], optional (table footnotes like [a], [b], etc.)
          - abbreviations: str, optional (abbreviation definitions)
    title : str, optional
        Optional header title for the notes file.

    Returns
    -------
    str
        Formatted notes text.
    """
    lines: List[str] = []

    if title:
        lines.append(title.upper())
        lines.append("=" * 78)
        lines.append(f"Generated: {_ts_iso()}")
        lines.append("")

    for i, legend in enumerate(legends):
        target = legend.get("target", f"Output {i + 1}")
        caption = legend.get("caption", "")
        footnotes = legend.get("footnotes", [])
        abbreviations = legend.get("abbreviations", "")

        # Target header
        lines.append(f"[{target}]")
        lines.append("")

        # Caption (wrapped)
        if caption:
            wrapped = textwrap.wrap(caption, width=78)
            lines.extend(wrapped)
            lines.append("")

        # Footnotes
        if footnotes:
            for fn in footnotes:
                fn_wrapped = textwrap.wrap(fn, width=74)
                lines.append(fn_wrapped[0])
                for cont in fn_wrapped[1:]:
                    lines.append("    " + cont)
            lines.append("")

        # Abbreviations
        if abbreviations:
            abbr_line = f"Abbreviations: {abbreviations}"
            abbr_wrapped = textwrap.wrap(abbr_line, width=78)
            lines.extend(abbr_wrapped)
            lines.append("")

        # Separator between entries
        lines.append("-" * 78)
        lines.append("")

    return "\n".join(line.rstrip() for line in lines).rstrip() + "\n"


# Helper formatters for dynamic values
def format_p(value: Optional[float], digits: int = 3) -> str:
    """Format p-value."""
    if value is None or not isinstance(value, (int, float)):
        return "p = NA"
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.{digits}f}"


def format_ci(estimate: float, lower: float, upper: float, digits: int = 2) -> str:
    """Format estimate with 95% CI."""
    return f"{estimate:.{digits}f} (95% CI: {lower:.{digits}f}–{upper:.{digits}f})"


def format_median_iqr(median: float, q1: float, q3: float, digits: int = 1) -> str:
    """Format median (IQR)."""
    return f"{median:.{digits}f} ({q1:.{digits}f}–{q3:.{digits}f})"


__all__ = [
    "build_notes",
    "format_p",
    "format_ci",
    "format_median_iqr",
]
