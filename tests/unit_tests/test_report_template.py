"""Focused tests for the shared report layout's web-only navbar additions."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from jinja2 import Environment, FileSystemLoader


TEMPLATES = Path(__file__).parents[2] / 'CRISPResso2' / 'CRISPRessoReports' / 'templates'


def test_layout_compiles_and_gates_plot_assistant_link():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    template = env.get_template('layout.html')
    source = template.environment.loader.get_source(env, 'layout.html')[0]

    assert 'current_user.is_authenticated and plot_assistant_enabled' in source
    assert 'plot_assistant_url or' in source
    assert 'AI Plot Assistant' in source


@pytest.mark.parametrize('child', [False, True])
@pytest.mark.parametrize('gate', ['enabled', 'disabled', 'default-user', 'anonymous', 'cli'])
def test_layout_assistant_targets_and_visibility(child, gate):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    context = {
        'is_web': gate != 'cli',
        'is_default_user': gate == 'default-user',
        'current_user': SimpleNamespace(is_authenticated=gate != 'anonymous', email='test@example.com', role='User'),
        'config': {'BANNER_TEXT': 'CRISPResso'},
        'get_flashed_messages': lambda **kwargs: [],
        'url_for': lambda endpoint: '/' + endpoint,
        'plot_assistant_enabled': gate != 'disabled',
        'plot_assistant_folder_id': 'parent',
    }
    expected = '/plot-assistant/parent'
    if child:
        expected += '/children/sample'
        context['plot_assistant_url'] = expected
    # Without an explicit URL, older C2Web versions still get a standalone link.
    html = env.get_template('layout.html').render(**context)
    if gate == 'enabled':
        assert f'href="{expected}"' in html
        assert html.count('AI Plot Assistant') == 1
        if child:
            assert 'href="/plot-assistant/parent"' not in html
    else:
        assert 'AI Plot Assistant' not in html


def test_core_report_header_uses_feature_flag_and_explicit_target():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)))
    source = env.loader.get_source(env, 'report.html')[0]
    assert "report_data.get('plot_assistant_enabled')" in source
    assert 'href="{{ report_data[\'plot_assistant_url\'] }}"' in source
    assert 'href="/plot-assistant/{{ report_data[\'folder_id\'] }}"' not in source


@pytest.mark.parametrize('default,authenticated', [(False, False), (True, True), (False, True)])
def test_web_chrome_does_not_change_with_identity(default, authenticated):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    html = env.get_template('layout.html').render(
        is_web=True, is_default_user=default, PUBLIC_BUILD=True,
        current_user=SimpleNamespace(is_authenticated=authenticated, email='user@example.edu', role='Admin'),
        config={'BANNER_TEXT': 'Legacy branding'},
        get_flashed_messages=lambda **kwargs: [], url_for=lambda endpoint: '/' + endpoint,
    )
    assert html.count('id="cr-site-nav"') == 1
    assert '<span>CRISPResso2</span>' in html
    assert 'class="crispresso-header"' not in html
    assert 'id="copyright_div"' in html
    signed_in = authenticated and not default
    assert ('href="/logout"' in html) == signed_in
    assert ('>Admin</a>' in html) == signed_in
    assert ('href="/login"' in html) == (not signed_in)


def test_standalone_chrome_has_no_web_dependencies():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    html = env.get_template('layout.html').render(is_web=False)
    assert 'class="crispresso-header"' in html
    assert 'id="cr-site-nav"' not in html
    assert 'Nat Biotechnol. 2019' in html
    assert 'Kendell Clement and Luca Pinello</a>' in html
