"""Focused tests for the shared report layout's web-only navbar additions."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from jinja2 import ChoiceLoader, DictLoader, Environment, FileSystemLoader


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
def test_web_chrome_brand_and_actions_follow_identity(default, authenticated):
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    html = env.get_template('layout.html').render(
        is_web=True, is_default_user=default, PUBLIC_BUILD=True,
        current_user=SimpleNamespace(is_authenticated=authenticated, email='user@example.edu', role='Admin'),
        config={'BANNER_TEXT': 'Legacy branding'},
        get_flashed_messages=lambda **kwargs: [], url_for=lambda endpoint: '/' + endpoint,
    )
    assert html.count('id="cr-site-nav"') == 1
    assert '<script src="/static/js/site_chrome.js" defer></script>' in html
    assert "<div id='left_help_div'>" in html
    signed_in = authenticated and not default
    brand = 'CRISPRessoPro' if signed_in else 'CRISPResso2'
    assert f'<span>{brand}</span>' in html
    assert 'class="crispresso-header"' not in html
    assert 'id="copyright_div"' in html
    assert ('href="/logout"' in html) == signed_in
    assert ('>Admin</a>' in html) == signed_in
    assert ('href="/login"' in html) == (not signed_in)


def test_standalone_chrome_has_no_web_dependencies():
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=True)
    html = env.get_template('layout.html').render(is_web=False)
    assert 'class="crispresso-header"' in html
    assert 'id="cr-site-nav"' not in html
    assert '/static/js/site_chrome.js' not in html
    assert 'Nat Biotechnol. 2019' in html
    assert 'Kendell Clement and Luca Pinello</a>' in html


LOGIN_PARTIAL = 'shared/partials/pro_login_button.html'
LOGIN_MARKER = '<div id="cr-pro-login-btn">STUB</div>'


def _layout_env(partials):
    return Environment(
        loader=ChoiceLoader([FileSystemLoader(str(TEMPLATES)), DictLoader(partials)]),
        autoescape=True,
    )


def _layout_context(is_web, is_default_user=False, authenticated=True):
    return {
        'is_web': is_web,
        'is_default_user': is_default_user,
        'current_user': SimpleNamespace(is_authenticated=authenticated, email='a@b.c', role='User'),
        'config': {'BANNER_TEXT': 'CRISPResso'},
        'get_flashed_messages': lambda **kwargs: [],
        'url_for': lambda endpoint: '/' + endpoint,
    }


@pytest.mark.parametrize('is_default_user, authenticated', [(True, True), (False, True), (False, False)])
def test_web_layout_ignores_obsolete_floating_login_partial(is_default_user, authenticated):
    html = _layout_env({LOGIN_PARTIAL: LOGIN_MARKER}).get_template('layout.html').render(
        **_layout_context(True, is_default_user, authenticated))
    assert 'cr-pro-login-btn' not in html
    assert html.count('id="cr-site-nav"') == 1


def test_web_layout_renders_without_login_partial():
    html = _layout_env({}).get_template('layout.html').render(**_layout_context(True, True))
    assert 'cr-pro-login-btn' not in html
    assert '</body>' in html


def test_cli_layout_never_includes_login_partial():
    with_partial = _layout_env({LOGIN_PARTIAL: LOGIN_MARKER}).get_template('layout.html').render(**_layout_context(False))
    without_partial = _layout_env({}).get_template('layout.html').render(**_layout_context(False))
    assert 'cr-pro-login-btn' not in with_partial
    assert with_partial == without_partial


def test_default_user_navbar_has_one_login_link():
    """An older C2Web floating partial must not duplicate the navbar sign-in."""
    button = '<a id="cr-pro-login-btn" href="/login">STUB</a>'
    html = _layout_env({LOGIN_PARTIAL: button}).get_template('layout.html').render(**_layout_context(True, True))
    assert html.count('href="/login"') == 1


@pytest.mark.parametrize('is_web', [False, True])
def test_layout_loads_twitter_widget_only_outside_the_web(is_web):
    # C2Web loads widgets.js itself, after cookie consent.
    html = _layout_env({}).get_template('layout.html').render(**_layout_context(is_web, is_default_user=True))
    assert ('twitter-share-button' in html) is (not is_web)
    assert ('platform.twitter.com/widgets.js' in html) is (not is_web)


def test_cli_layout_twitter_loader_bytes_are_unchanged():
    # CRISPResso2_tests diffs report HTML against baselines, so the web gate must not move a byte here.
    import logging

    from CRISPResso2.CRISPRessoReports.CRISPRessoReport import get_jinja_loader

    env = get_jinja_loader(str(TEMPLATES.parents[1]), logging.getLogger(__name__))
    html = env.get_template('layout.html').render(is_web=False, is_default_user=False)
    assert '</script>\n\n\n  <script>!function(d,s,id){' in html
    assert "(document, 'script', 'twitter-wjs');</script>\n  \n\n  </head>" in html


def test_layout_loads_twitter_widget_when_is_web_is_not_passed():
    html = _layout_env({}).get_template('layout.html').render()
    assert html.count('platform.twitter.com/widgets.js') == 1
