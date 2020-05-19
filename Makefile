ATOMPACK_DIR=$(shell pwd)
MYPY_DIR=$(ATOMPACK_DIR)/mypy

clean:
	@find . | grep -E "(__pycache__|\.pyc)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf

format:
	@yapf --in-place --recursive --parallel --style="google" $(ATOMPACK_DIR)

lint:
	@export MYPYPATH=$(MYPY_DIR);\
		mypy --config-file=$(MYPY_DIR)/mypy.ini $(ATOMPACK_DIR)/atompack/
	@pyflakes $(ATOMPACK_DIR)/atompack

test:
	@pytest --doctest-modules -v $(ATOMPACK_DIR)/atompack/
	@make clean
