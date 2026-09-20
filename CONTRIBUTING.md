# Contributing to swesphere

Bug reports, questions and pull requests are welcome.

## Reporting a problem

Open an issue at https://github.com/enino84/swesphere/issues with:

- the version (`python -c "import swesphere; print(swesphere.__version__)"`), your Python version and operating system;
- a short script that reproduces the problem, and what you expected instead.

For questions about using the model, the same issue tracker is fine, or write to enino@uninorte.edu.co.

## Proposing a change

1. Fork the repository and create a branch.
2. Install in editable mode with the test extras: `pip install -e ".[test]"`.
3. Make the change and add a test for it in `tests/`.
4. Run `python -m pytest tests -q`. All tests must pass.
5. Add a line to `CHANGELOG.md` and open a pull request describing what changed and why.

Changes to the numerics (the right-hand side, the filter, the sponge, the forcing or the shipped
integrator) alter the statistics of the documented regimes. If your change touches them, rerun the
relevant experiment (`docker compose run --rm verification`, `stability` or `regimes`) and report
the numbers in the pull request.

## Adding a time-stepping scheme

A scheme is a function `step(rhs, u, v, h, dt) -> (u, v, h)` registered in
`swesphere.integrators.INTEGRATORS`. See the README for an example. A scheme proposed for inclusion
should come with its stable range of `dt` from EXP-01.

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
