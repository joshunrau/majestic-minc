#AUTOMAKE_OPTIONS = check-news

EXTRA_DIST = doxymain.txt

include_HEADERS = arguments.h

LDADD = libarguments.la

lib_LTLIBRARIES = libarguments.la

TESTS = test

check_PROGRAMS = $(TESTS)

libarguments_la_SOURCES = arguments.cpp

# Tests
test_SOURCES = test.cpp






