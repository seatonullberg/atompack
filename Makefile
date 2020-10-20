MYPYPATH=./mypy

bench:
	@pipenv run pytest --benchmark-columns="min, median, max, stddev, rounds, iterations" -v ./benches 

clean:
	@find . | grep -E "(.benchmarks)" | xargs rm -rf
	@find . | grep -E "(.cache)" | xargs rm -rf
	@find . | grep -E "(.mypy_cache)" | xargs rm -rf
	@find . | grep -E "(.pytest_cache)" | xargs rm -rf
	@find . | grep -E "(__pycache__)" | xargs rm -rf
	@find . | grep -E "(build)" | xargs rm -rf
	@find . | grep -E "(dist)" | xargs rm -rf
	@find . | grep -E "(\.pyc|\.so|\.egg-info)" | xargs rm -rf

document:
	@pipenv run pdoc --html --force\
		--template-dir ./docs/config\
		--output-dir ./docs\
		./atompack

format:
	@pipenv run isort ./atompack
	@pipenv run yapf -rip --style='{based_on_style: google, column_limit: 120}' ./

lint:
	@export MYPYPATH=$(MYPYPATH);\
		pipenv run mypy --config-file=$(MYPYPATH)/mypy.ini ./atompack
	@pipenv run pyflakes ./atompack

test:
	@pipenv run pytest --doctest-modules -v ./atompack ./tests
