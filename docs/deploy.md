# Streamlit Community Cloud deployment

This guide prepares the app for deployment. The app has not been deployed from
this checkout.

## Deployment settings

- **GitHub repository:** Unverified. This checkout has no configured Git
  remote. Select the GitHub repository containing this project in Community
  Cloud.
- **Branch:** Unverified on GitHub. The current local branch is `main`; select
  the matching branch only after confirming it exists in the chosen repository.
- **Main file path:** `src/colony_breeder/app/streamlit_app.py`
- **Requirements file:** root `requirements.txt`, containing `-e .[app]`.
  This installs this project as an editable package with its `[app]` extra.
  That extra includes Streamlit, pandas, openpyxl, and defusedxml. The
  setuptools configuration discovers the `colony_breeder` package under
  `src/`, so the page can import it after package installation.
- **Python version recommendation:** Python 3.12. This is a recommendation,
  not a change to the project's `requires-python = ">=3.9"` setting. The
  current Streamlit Community Cloud documentation says its default is Python
  3.12. The app has not been tested in this checkout on Python 3.12; local
  tests use Python 3.14.

## Deploy using Community Cloud

1. Confirm the intended repository and branch contain the app, the root
   `requirements.txt`, `.streamlit/config.toml`, and only fake demonstration
   data.
2. Open your Streamlit Community Cloud workspace and choose **Create app**.
3. Select the confirmed repository and branch, then enter
   `src/colony_breeder/app/streamlit_app.py` as the main file path.
4. In **Advanced settings**, select the recommended Python 3.12 version.
5. Review the settings and choose **Deploy**. Watch the build and app logs for
   installation or startup errors.

The repository and branch names are unverified because no remote is configured
in this checkout. No live Cloud build or deployment was attempted, so runtime
installation and startup on Community Cloud remain unverified.

## Upload data warning

A hosted app processes uploaded files on the host's server. A public repository
or public app must never contain real lab data. Do not commit real lab
spreadsheets or use real colony data in a public deployment. The app's upload
note describes its in-memory handling; it does not mean hosted uploads are
processed on the user's computer.
