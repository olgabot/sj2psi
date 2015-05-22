test:
	py.test

coverage:
	coverage run --source sj2psi --omit=test --module py.test

lint:
	pyflakes sj2psi

pep8:
	pep8 sj2psi