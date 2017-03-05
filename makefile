go:
	go build
	-./go-consensus -h
	${MAKE} test
clean:
	go clean
uninstall:
	rm -f ../../bin/go-consensus
test:
	make -C t/ new
test-old:
	make -C t/ old
