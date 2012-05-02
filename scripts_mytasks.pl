#!/usr/bin/perl

use strict;

my $date = `date`;
print "# $date";

sub escape {
    my $text = shift;
    $text =~ s/\\/\\\\/g;
    $text =~ s/\n/\\n/g;
    $text =~ s/\t/\\t/g;
    return $text;
}
sub checkLine {
    my $line = shift;
    my $file = shift;
    my $row = shift;

    $file = escape($file);
    $row = escape($row);

    if ($file =~ /(.cpp|.c|.cxx|.h)$/ && $line =~ /(FIXME|TODO)(.*)$/) {
        print "$file\t$row\tTODO\t$1$2\n";
    }
    if ($line =~ /^(<<<<|====|>>>>)/) {
        print "$file\n$row\tERROR\tResolved conflict.\n";
    }
    if ($file =~ /(.cpp|.c|.cxx|.h)$/ && $line =~ /(Q_SIGNALS|Q_SLOTS)/) {
        print "$file\t$row\tWARNING\tPrefer signals and slots over Q_SIGNALS and Q_SLOTS.\n";
    }
    if ($line =~ /^#\s*(warning|error) (.*)$/) {
        print "$file\t$row\t$1\tClean up preprocessor $1:$2.\n";
    }
    if ($file =~ /(.cpp|.c|.cxx|.h)$/ && $line =~ /\s+$/) {
#        print "$file\t$row\tWARNING\tTrailing whitespaces found.\n";
    }
    if ($file =~ /(.cpp|.c|.cxx|.h)$/ && $line =~ /^ *\t */) {
        print "$file\t$row\tWARNING\tTabs used to indent.\n";
    }
    if ($file =~ /(.cpp|.c|.cxx|.h)$/ && $line =~ /{.*[^\s].*$/) {
        if ($line !~ /{\s+\/\*/)
        { print "$file\t$row\tWARNING\tText found after opening bracket.\n"; }
    }
}

my $filename = "";
my $pos = 0;

sub getDiffOrigin {
   my $currentBranch = `git branch | grep "^* "`;
   chop $currentBranch;
   $currentBranch =~ s/^* //;

   return "958504f9eedf0d9c8ff65435dab4dad4d83e1fa0";

   my $remoteBranch = `git config --local --get branch.$currentBranch.merge`;
   chop $remoteBranch;
   $remoteBranch =~ s!^refts/heads/!!;

   return "958504f9eedf0d9c8ff65435dab4dad4d83e1fa0" if (!$remoteBranch);
}

my $origin = getDiffOrigin;
print "# running: git diff $origin ...\n";

open(PIPE, "git diff $origin|");
while (<PIPE>) {
    chomp;
    my $line = $_;
    print "#$line\n";
    next if ($line =~ /^-/);
    if ($line =~ /^\+\+\+ (.*)$/) {
        $filename = $1;
        $filename =~ s/^b\///;
        $pos = 0;
        next;
    }
    next if $filename =~ /^\/dev\/null$/;
    if ($line =~ /^@@ -\d+,\d+\s\+(\d+),\d+ @@$/) {
        $pos = $1 - 1;
        next;
    }
    $pos = $pos + 1;
    if ($line =~ /^\+(.*)/) {
        checkLine($1, $filename, $pos);
    }
}  
close(PIPE);  
