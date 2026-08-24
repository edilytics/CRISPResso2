"""Render-contract tests for the shared hover-zoom figure macro.

``shared/partials/figure_zoom.html`` is the single hover-zoom implementation
for report figures (legacy plot_2a): it is rendered by the no-Pro
CRISPRessoReports report template AND by CRISPRessoPro's nucleotide-quilt
partial, so its markup and inline behavior are pinned here. The bugs that
motivated the shared implementation were all markup/JS regressions:

- the zoom strip collapsing to 0px wide inside flex-column wrappers
  (missing explicit width)
- one-shot JS state ("hasWidthSet") locking in zero/NaN geometry when the
  first hover happened in a hidden tab or before image load
- Pro and no-Pro renderers drifting apart (different zoom markup, zoom on
  figures where the zoomed copy renders smaller than the plain figure)

The macro must render identically under every jinja environment the report
stack uses: the CLI builder and CRISPRessoPro's hooks (autoescape off) and
the CRISPRessoWEB Flask app (autoescape by file extension). The inline
script is included from figure_zoom.js and must stay unescaped everywhere —
a report zip opened from file:// cannot load external scripts.
"""

import os

import pytest
from jinja2 import Environment, FileSystemLoader, select_autoescape

from CRISPResso2 import CRISPRessoReports

TEMPLATES_DIR = os.path.join(os.path.dirname(CRISPRessoReports.__file__), "templates")
JS_PATH = os.path.join(TEMPLATES_DIR, "shared", "partials", "figure_zoom.js")


def _env(flask_like):
    loader = FileSystemLoader(TEMPLATES_DIR)
    if flask_like:
        return Environment(loader=loader, autoescape=select_autoescape(enabled_extensions=("html", "htm", "xml")))
    return Environment(loader=loader)


def _render(flask_like=False, **macro_kwargs):
    kwargs = {"img_src": "out/2a.REF1.Nucleotide_percentage_quilt.png",
              "pdf_href": "out/2a.REF1.Nucleotide_percentage_quilt.pdf",
              "uid": "nucs_REF1"}
    kwargs.update(macro_kwargs)
    env = _env(flask_like)
    template = env.from_string(
        "{% from 'shared/partials/figure_zoom.html' import figure_zoom %}"
        "{{ figure_zoom(img_src, pdf_href, uid, width) }}"
    )
    return template.render(width="95%", **kwargs)


@pytest.mark.parametrize("flask_like", [False, True], ids=["cli/pro-env", "flask-env"])
class TestFigureZoomMarkup:
    def test_strip_has_explicit_width_and_overflow_window(self, flask_like):
        """The strip must carry width:100% + overflow:hidden.

        It is a flex item inside flex-column figure wrappers; without an
        explicit width it shrink-to-fits to 0px (the blank-strip bug) and
        without overflow:hidden the translated image spills instead of
        cropping.
        """
        html = _render(flask_like)
        assert 'data-cr-zoom-strip="nucs_REF1"' in html
        strip_style = html.split('data-cr-zoom-strip="nucs_REF1"', 1)[1].split(">", 1)[0]
        assert "width:100%" in strip_style.replace(" ", "")
        assert "overflow:hidden" in strip_style.replace(" ", "")

    def test_all_five_elements_carry_matching_data_attributes(self, flask_like):
        """strip / strip-img / figure / lens / img are bound by data-cr-zoom-*."""
        html = _render(flask_like)
        for attr in ("data-cr-zoom-strip=", "data-cr-zoom-strip-img=",
                     "data-cr-zoom-figure=", "data-cr-zoom-lens=", "data-cr-zoom-img="):
            assert attr in html, f"missing {attr} binding"

    def test_pdf_link_wraps_main_image(self, flask_like):
        html = _render(flask_like)
        assert '<a href="out/2a.REF1.Nucleotide_percentage_quilt.pdf">' in html
        assert 'src="out/2a.REF1.Nucleotide_percentage_quilt.png"' in html

    def test_inline_script_is_the_shared_stateless_js(self, flask_like):
        """The figure embeds figure_zoom.js verbatim (reports are self-contained).

        Unescaped in every environment: escaping would break require()-based
        unit tests parity and, more importantly, the browser script itself.
        """
        html = _render(flask_like)
        js = open(JS_PATH).read()
        assert js in html, "figure_zoom.js must be inlined verbatim by the macro"

    def test_script_is_stateless_no_cached_geometry(self, flask_like):
        """The old one-shot init ('hasWidthSet') bug must not come back."""
        html = _render(flask_like)
        assert "lens.hasWidthSet" not in html
        assert "function computeZoom" in html
        assert "pointermove" in html

    def test_strip_images_defy_bootstrap_img_clamping(self, flask_like):
        """Bootstrap's img{max-width:100%} would clamp and distort the
        full-height strip copy; the macro must opt out explicitly."""
        html = _render(flask_like)
        assert html.count("max-width:none") >= 2  # zoom copy + tablet scroll copy

    def test_tablet_scroll_fallback_preserved(self, flask_like):
        """Below lg there is no hover: the strip doubles as a scrollable view."""
        html = _render(flask_like)
        assert 'class="d-lg-none"' in html
        assert "overflow-x:auto" in html.replace(" ", "")


def test_uid_is_escaped_in_flask_env():
    """A hostile uid must not break out of the attribute in the Flask env.

    (The CLI/Pro jinja environments run without autoescape — same trust model
    as the pre-existing report code: fig roots come from run output, not user
    input. The internet-facing path is C2Web's Flask env, which escapes.)
    """
    html = _render(flask_like=True, uid='a"><script>')
    escaped = 'a&#34;&gt;&lt;script&gt;'
    assert escaped in html, 'uid must be HTML-escaped in the Flask env'
    assert html.count('data-cr-zoom-strip="' + escaped) == 1
    assert html.count('data-cr-zoom-figure="' + escaped) == 1



class TestNucQuiltPartial:
    """``partials/nuc_quilt_zoom_figure.html`` is the quilt shim around the
    shared macro (the drop-in for CRISPRessoPro's partial). Its four
    rendering branches and the zoom-precedence rule are pinned here because
    a wrong branch silently changes the report figure (e.g. zooming an
    sgRNA-scoped quilt renders it SMALLER than the plain figure)."""

    def _render_fig(self, **overrides):
        fig = {
            "id": "2a", "root": "2a.REF1.Nucleotide_percentage_quilt",
            "amplicon": "REF1", "sgRNA_index": None, "caption": "cap",
            "data_files": [("txt", "f.txt")], "width": "100%",
        }
        fig.update(overrides)
        template = _env(False).from_string(
            "{% include 'partials/nuc_quilt_zoom_figure.html' %}"
        )
        return template.render(fig=fig, data_path="out/")

    def test_json_data_branch_renders_d3_widget(self):
        html = self._render_fig(json_data='{"a":1}')
        assert 'id="nuc_quilt_2a_REF1_Nucleotide_percentage_quilt"' in html
        assert 'const nuc_quilt_2a_REF1_Nucleotide_percentage_quilt = {"a":1};' in html
        assert "data-cr-zoom-" not in html

    def test_html_branch_renders_fragment_verbatim(self):
        html = self._render_fig(html="<b>pre</b>")
        assert "<b>pre</b>" in html
        assert "data-cr-zoom-" not in html

    def test_amplicon_wide_quilt_without_zoom_key_zooms(self):
        html = self._render_fig()  # sgRNA_index None, no zoom key -> fallback
        assert 'data-cr-zoom-figure="2a_REF1_Nucleotide_percentage_quilt"' in html
        assert 'src="out/2a.REF1.Nucleotide_percentage_quilt.png"' in html

    def test_sgrna_scoped_quilt_renders_plain_image(self):
        html = self._render_fig(sgRNA_index=1)
        assert "width='100%'" in html  # plain image + pdf link, no zoom strip
        assert "data-cr-zoom-" not in html
        assert 'id="fig_2a_REF1_sg1"' in html  # sgRNA suffix on the wrapper id

    def test_explicit_fig_zoom_overrides_scope_fallback(self):
        assert "data-cr-zoom-" not in self._render_fig(zoom=False)
        assert "data-cr-zoom-figure" in self._render_fig(sgRNA_index=0, zoom=True)

    def test_missing_optional_keys_do_not_crash(self):
        """Older CRISPRessoPro figure dicts may omit id/caption/data_files.

        Iterating an undefined ``data_files`` raises UndefinedError even in
        non-strict jinja environments, which would kill the whole report for
        one figure. A minimal fig (root only) must render, not raise.
        """
        template = _env(False).from_string(
            "{% include 'partials/nuc_quilt_zoom_figure.html' %}"
        )
        html = template.render(
            fig={"root": "2a.REF1.Nucleotide_percentage_quilt"}, data_path="out/"
        )
        assert 'id="fig_"' in html  # missing id defaults to empty, not an error
        assert "Data:" not in html  # missing data_files renders no data links
        assert 'src="out/2a.REF1.Nucleotide_percentage_quilt.png"' in html


class TestReportZoomWiring:
    """report.html must wire plot_2a locs through the shared macro (and must
    not resurrect the removed per-amplicon zoom ids/JS it used to inline)."""

    def test_plot_2a_locs_render_shared_zoom(self):
        from CRISPResso2.CRISPRessoReports.jinja_partials import (
            generate_render_partial, render_partial,
        )

        env = _env(False)

        def custom_partial_render(name, **data):
            template = env.get_template(name)
            data.update(
                render_partial=generate_render_partial(custom_partial_render),
                is_default_user=False, is_web=False, C2PRO_INSTALLED=False,
            )
            return template.render(**data)

        report_data = {
            "amplicons": ["REF1"],
            "figures": {
                "htmls": {"REF1": {}},
                "locs": {"REF1": {"plot_2a": "2a.REF1.Nucleotide_percentage_quilt"}},
                "titles": {"REF1": {"plot_2a": "Figure 2a"}},
                "captions": {"REF1": {"plot_2a": "cap"}},
                "datas": {"REF1": {"plot_2a": []}},
                "names": {"REF1": ["plot_2a"]},
                "sgRNA_based_names": {"REF1": {}},
            },
            "run_data": {"running_info": {"args": {}},
                         "results": {"refs": {"REF1": {}}, "ref_names": ["REF1"],
                                     "general_plots": {}}},
            "report_display_name": "test", "crispresso_data_path": "",
            "nuc_quilt_names": [], "allele_table_names": [],
        }
        html = render_partial("report.html", custom_partial_render,
                              report_data=report_data, C2PRO_INSTALLED=False)
        assert 'data-cr-zoom-figure="nucs_REF1"' in html
        assert 'src="2a.REF1.Nucleotide_percentage_quilt.png"' in html
        assert 'href="2a.REF1.Nucleotide_percentage_quilt.pdf"' in html
        assert html.count("function computeZoom") == 1
        # the retired per-amplicon zoom machinery must stay gone
        for gone in ("tozoom_nucs_REF1", "zoomview_nucs_REF1",
                     "zoomlens_nucs_REF1", "updateZoom"):
            assert gone not in html
