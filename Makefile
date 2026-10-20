ENV_NAME := rosetta_pipelines

.PHONY: help env install-pyrosetta install test format clean

help:
	@echo "Targets:"
	@echo "  make env               Create the conda env, install PyRosetta, and install this package"
	@echo "  make install-pyrosetta Install PyRosetta into the active env (via pyrosetta-installer)"
	@echo "  make install           pip install this package into the active env"
	@echo "  make test              Run the test suite"
	@echo "  make format            Format with black (line length 99)"
	@echo "  make clean             Remove build artifacts"

# One-shot setup: build the env from environment.yml (which pip-installs pyrosetta-installer and this
# package), then run the installer to fetch the actual PyRosetta wheel for this platform/Python.
env:
	conda env create -f environment.yml
	conda run -n $(ENV_NAME) python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
	@echo "Done. Activate with:  conda activate $(ENV_NAME)"

# Install PyRosetta into whatever env is currently active (RosettaCommons official helper).
install-pyrosetta:
	pip install pyrosetta-installer
	python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

# Install this package into the active env.
install:
	pip install .

test:
	pytest tests/

format:
	black --line-length 99 rosetta_pipelines tests

clean:
	rm -rf build dist *.egg-info
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
