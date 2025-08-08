# Contributing guidelines for _Campylobacter_ CRISPRscape

Focused on development of analysis code and documentation.

## Main workflow

We aim to follow a method of working that:

1. Allows different people to work on separate features without breaking
each other's code

2. Moves development away from the main branch, preventing users from
cloning a non-functional workflow

3. Has a simple release model (with version numbers and release notes)

Therefore, starting 8 August 2025, the general git-based working procedure is as follows:

- Next to the main branch, there is a 'dev' branch.

- From the dev branch, each contributor makes a new branch with the
name of the feature they're working on. (For example, 'cctyper'.)

- The person works on this branch until the feature appears to work.

- The function-focused branch is merged into dev and pushed to GitHub.

- The dev branch is pulled to a test folder with small dataset
(this is currently done by @samnooij).

- Only when the Snakemake workflow finished completely and successfully,
the dev branch is merged into the main branch and pushed to GitHub.

- Once a new version of the main branch is pushed, a tagged release
is created with a short description of the solved problem/new feature.

## Language-specific practices

This project uses a combination of programming and markup languages.
For most of them, we try to follow styling conventions.
These include:

- [Black](https://github.com/psf/black) for Python

- [snakefmt](https://github.com/snakemake/snakefmt) for Snakemake

- [styler](https://styler.r-lib.org/) for R

- [markdownlint](https://github.com/DavidAnson/markdownlint/) for Markdown

**Please format code and Markdown files before committing to git.**
(Some editors, such as [vscodium](https://vscodium.com/) have extensions that
can automatically apply linting rules upon saving!)

## Documentation

Documentation is written both in the `README.md` file, and in separate
Markdown files under `doc/`. The former functions as GitHub homepage and
should at least cover the basics about this project. More detailed information
is given in the separate Markdown files under `doc/`, which are automatically
converted to nice GitHub pages. (Provided that they are correctly linked in
the file `mkdocs.yaml`.))
