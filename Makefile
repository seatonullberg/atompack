MYPYPATH=./mypy

bench:
	@make build
	@python3 -m pytest -v ./atompack ./benches 

build:
	@python3 setup.py build_ext --inplace

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
	@python3 -m isort -rc ./atompack
	@python3 -m yapf -rip --style='{based_on_style: google, column_limit: 120}' ./

lint:
	@export MYPYPATH=$(MYPYPATH);\
		python3 -m mypy --config-file=$(MYPYPATH)/mypy.ini ./atompack
	@python3 -m pyflakes ./atompack

test:
	@make build
	@python3 -m pytest --doctest-modules -v ./atompack ./tests
