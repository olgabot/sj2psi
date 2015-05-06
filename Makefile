test:
	py.test

coverage:
	py.test --cov sj2psi --cov-report term-missing sj2psi/test

lint:
	pyflakes sj2psi

pep8:
	pep8 sj2psi