"""Batch and Aggregate share modification tabs and data-link behavior."""

import logging
import re
from pathlib import Path

import pytest

from CRISPResso2.CRISPRessoReports import CRISPRessoReport


@pytest.mark.parametrize('mode', ['batch', 'aggregate'])
@pytest.mark.parametrize('plot_count', [0, 1, 8])
@pytest.mark.parametrize('report_inside', [False, True])
def test_modification_summary_tabs(tmp_path, monkeypatch, mode, plot_count, report_inside):
    monkeypatch.setattr(CRISPRessoReport, 'C2PRO_INSTALLED', False)
    output = tmp_path / 'output'
    output.mkdir()
    report = (output if report_inside else tmp_path) / 'report.html'
    plots = {}
    for kind in ('heatmap', 'line'):
        names = [f'{kind}-{i}' for i in range(plot_count)]
        plots[kind] = {
            'names': names,
            'htmls': {name: f'<div id="plot-{name}" class="plotly-graph-div"></div>' for name in names},
            'titles': {name: f'Title {name}' for name in names},
            'labels': {name: f'Caption {name}' for name in names},
            'datas': {name: [('Modification fractions', f'{name}.txt')] for name in names},
        }
    CRISPRessoReport.make_multi_report(
        ['child'], [], [], {'child': 'child.html'}, str(report), str(output),
        str(Path(CRISPRessoReport.__file__).parents[1]), 'Test report', mode,
        logging.getLogger(__name__),
        allele_modification_heatmap_plot=plots['heatmap'],
        allele_modification_line_plot=plots['line'],
    )
    html = report.read_text()
    prefix = '' if report_inside else 'output/'
    assert f'href="{prefix}child.html"' in html
    assert ('<h5>Summary Plots</h5>' in html) == bool(plot_count)
    ids = re.findall(r'\bid="([^"]+)"', html)
    assert len(ids) == len(set(ids))
    targets = re.findall(r'data-bs-target="#([^"]+)"', html)
    assert len(targets) == plot_count * 2
    for i in range(plot_count):
        for kind in ('heatmap', 'line'):
            pane = f'{mode}-modification-{i}-{kind}'
            assert pane in targets
            assert pane in ids
            assert f'aria-controls="{pane}"' in html
            assert f'aria-labelledby="{pane}-tab"' in html
            assert f'id="plot-{kind}-{i}"' in html
            assert f'Caption {kind}-{i}' in html
            assert f'href="{prefix}{kind}-{i}.txt"' in html
    if plot_count or mode == 'aggregate':
        assert f"#{mode}-modification-summary" in html
        assert "tab.addEventListener('shown.bs.tab'" in html
        assert "pane.querySelectorAll('.plotly-graph-div')" in html
        assert 'Plotly.Plots.resize(plot)' in html
