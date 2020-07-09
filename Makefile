MYPYPATH=./mypy

clean:
	@find . | grep -E "(.benchmarks)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(__pycache__)" | xargs rm -rf
	@find . | grep -E "(build)" | xargs rm -rf
	@find . | grep -E "(dist)" | xargs rm -rf
	@find . | grep -E "(\.pyc|\.so)" | xargs rm -rf

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
	@export MYPYPATH=$(MYPYPATH);\
		mypy --config-file=$(MYPYPATH)/mypy.ini ./atompack
	@pyflakes ./atompack

test:
	@python3 setup.py build_ext --inplace
	@pytest --doctest-modules -v ./atompack ./tests
