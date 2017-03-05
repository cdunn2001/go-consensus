go:
	go build
	./go-consensus
	${MAKE} test
clean:
	go clean
uninstall:
	rm -f ../../bin/go-consensus
test:
	make -C t/ new
test-old:
	make -C t/ old
