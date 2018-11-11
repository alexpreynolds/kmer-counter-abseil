SHELL=/bin/bash

#all: generate_fasta generate_test_data

all: kfs msl generate_fasta generate_test_data

abseil-cpp:
	git clone https://github.com/abseil/abseil-cpp.git

abseil-cpp-build: abseil-cpp
	cd abseil-cpp && bazel test //absl/...

generate_test_data: bazel-bin/src/generate_fasta
	bazel-bin/src/generate_fasta > test.fa

kfs: abseil-cpp
	cd src && bazel build //src:kfs --verbose_failures

msl: abseil-cpp
	cd src && bazel build //src:msl --verbose_failures

generate_fasta: abseil-cpp
	cd src && bazel build //src:generate_fasta --verbose_failures

clean:
	rm -f *~
	rm -f generate_fasta
	rm -f kfs
	rm -f msl
	rm -f test.fa