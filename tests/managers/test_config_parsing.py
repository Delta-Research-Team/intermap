"""
Creator: rglez
Date: 3/22/26
Env:
Desc:
"""

import intermap.managers.config as cfg_mod
from intermap.managers.config import print_colored_ascii


class TestASCIILogo:
    def test_execution(self, capsys):
        print_colored_ascii()
        assert "Error when reading HTML" not in capsys.readouterr().out

    def test_file_error(self, capsys):
        orig = cfg_mod.proj_dir
        cfg_mod.proj_dir = "/non/existent/path"
        try:
            print_colored_ascii()
            assert "Error when reading HTML" in capsys.readouterr().out
        finally:
            cfg_mod.proj_dir = orig
