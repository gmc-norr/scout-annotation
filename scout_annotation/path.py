import os
from pathlib import Path
import pathlib
import re


class WildcardPath:
    _wildcard_pattern: re.Pattern = re.compile(r"{([^}]+)}")
    _wildcard_sub_pattern: str = r"(?P<\1>.*?)"
    _forbidden_wildcard_chars: str = r"/\\"

    def __init__(self, path: str | Path):
        self.path: Path = Path(path)
        self.wildcard_matches: list[tuple[str, int, int]] = self._extract_wildcard_keys()
        if not self.wildcards:
            raise ValueError("no wildcards in path, use pathlib.Path instead")

    def __repr__(self) -> str:
        return f"<WildcardPath {str(self)}>"

    def __str__(self) -> str:
        return str(self.path)

    def _extract_wildcard_keys(self) -> list[tuple[str, int, int]]:
        return [
            (m.group(1), m.start(), m.end())
            for m in self._wildcard_pattern.finditer(str(self.path))
        ]

    @property
    def wildcards(self):
        if not self.wildcard_matches:
            return []
        return sorted(list(set([x[0] for x in self.wildcard_matches])))

    def determined_root(self):
        first_wildcard_start = self.wildcard_matches[0][1]
        prefix = str(self.path)[:first_wildcard_start]
        if not prefix.endswith(os.sep):
            return Path(prefix).parent
        return Path(prefix)

    def expand(self) -> list[tuple[Path, dict[str, str]]]:
        paths = []

        path = str(self.path)
        regex_parts = []
        var_positions = []
        last_end = 0

        for name, start, end in self.wildcard_matches:
            if start > last_end:
                regex_parts.append(re.escape(path[last_end:start]))
            regex_parts.append(rf"([^{self._forbidden_wildcard_chars}]*?)")
            var_positions.append(name)
            last_end = end

        if last_end < len(path):
            regex_parts.append(path[last_end:])

        extract_pattern = re.compile("".join(regex_parts))

        for root, dirnames, filenames in self.determined_root().walk():
            for fn in filenames:
                valid_wildcards = True
                wildcards = {}
                match = extract_pattern.fullmatch(str(root / fn))
                if match is None:
                    continue
                for i, name in enumerate(var_positions):
                    value = match.group(i + 1)
                    if name not in wildcards:
                        wildcards[name] = value
                    elif wildcards[name] != value:
                        valid_wildcards = False
                if valid_wildcards:
                    paths.append((root / fn, wildcards))

        return paths
