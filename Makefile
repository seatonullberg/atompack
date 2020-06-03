PWD=$(shell pwd)
MYPY_DIR=$(PWD)/mypy

clean:
	@find . | grep -E "(__pycache__|\.pyc)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf

document:
	@pdoc --html --force\
		--template-dir $(PWD)/docs\
		--output-dir $(PWD)/docs\
		$(PWD)/atompack
	@mv $(PWD)/docs/atompack/* $(PWD)/docs
	@rm -rf $(PWD)/docs/atompack

format:
	@isort $(PWD)/atompack/*.py
	@yapf -rip --style='{based_on_style: google, column_limit: 120}' $(PWD)/atompack

lint:
	@export MYPYPATH=$(MYPY_DIR);\
		mypy --config-file=$(MYPY_DIR)/mypy.ini $(PWD)/atompack/
	@pyflakes $(PWD)/atompack

publish:
	@python3 setup.py sdist bdist_wheel
	@twine check dist/*
	@twine upload dist/*

test:
	@pytest --doctest-modules -v $(PWD)/atompack/
	@make clean
