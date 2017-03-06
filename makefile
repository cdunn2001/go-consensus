go:
	-./go-consensus -h
	${MAKE} test
cgo:
	go tool cgo main.go
build:
	go build
clean:
	go clean
uninstall:
	rm -f ../../bin/go-consensus
test:
	make -C t/ new
test-old:
	make -C t/ old
