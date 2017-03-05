go:
	go build
	./go-falcon
	${MAKE} test
clean:
	go clean
uninstall:
	rm -f ../../bin/go-falcon
test:
	make -C t/ new
test-old:
	make -C t/ old
