.PHONY: reproduce test

reproduce:
	python reproduce.py

test:
	python -m unittest discover -s experiments -p 'test_*.py' -v
