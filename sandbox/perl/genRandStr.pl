#!/usr/bin/env perl
@char = ("A".."Z");
for (1..100_000) {
    for (0..10) {
	print $char[rand 5];
    }
    for (0..int(rand 90)) {
	print $char[rand @char];
    }
    print "\n";
}
