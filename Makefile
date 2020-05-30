ATOMPACK_DIR=$(shell pwd)
MYPY_DIR=$(ATOMPACK_DIR)/mypy

clean:
	@find . | grep -E "(__pycache__|\.pyc)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf

document:
	@pdoc --html --force\
		--template-dir $(ATOMPACK_DIR)/docs\
		--output-dir $(ATOMPACK_DIR)/docs\
		$(ATOMPACK_DIR)/atompack
	@mv $(ATOMPACK_DIR)/docs/atompack/* $(ATOMPACK_DIR)/docs
	@rm -rf $(ATOMPACK_DIR)/docs/atompack

format:
	@isort $(ATOMPACK_DIR)/atompack/*.py
	@yapf -rip --style='{based_on_style: google, column_limit: 120}' $(ATOMPACK_DIR)

lint:
	@export MYPYPATH=$(MYPY_DIR);\
		mypy --config-file=$(MYPY_DIR)/mypy.ini $(ATOMPACK_DIR)/atompack/
	@pyflakes $(ATOMPACK_DIR)/atompack

test:
	@pytest --doctest-modules -v $(ATOMPACK_DIR)/atompack/
	@make clean
