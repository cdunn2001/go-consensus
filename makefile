go:
	go build
	./go-falcon
clean:
	go clean
uninstall:
	rm -f ../../bin/go-falcon
test:
	make -C t/ old
