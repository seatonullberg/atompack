MYPY_DIR=./mypy

# TODO: build libatompack.so here
build:
	@python3 setup.py sdist bdist_wheel

clean:
	@find . | grep -E "(__pycache__|\.pyc)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf
	@find . | grep -E "(.benchmarks)" | xargs rm -rf

document:
	@pdoc --html --force\
		--template-dir ./docs\
		--output-dir ./docs\
		./atompack
	@mv ./docs/atompack/* ./docs
	@rm -rf ./docs/atompack

format:
	@isort -rc ./atompack
	@yapf -rip --style='{based_on_style: google, column_limit: 120}' ./atompack

lint:
	@export MYPYPATH=$(MYPY_DIR);\
		mypy --config-file=$(MYPY_DIR)/mypy.ini ./atompack
	@pyflakes ./atompack

publish:
	@make build
	@twine check dist/*
	@twine upload dist/*

test:
	@pytest --doctest-modules -v ./atompack
