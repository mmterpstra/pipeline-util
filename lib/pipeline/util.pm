package pipeline::util;

use 5.020002;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use pipeline-util ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = "0.8.5";


# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

pipeline::util - Perl cmdline tools for hts pipeline

=head1 SYNOPSIS

execute scripts from bash

=head1 DESCRIPTION

read readme on github

=head2 EXPORT

None by default.

=head1 SEE ALSO

	https://github.com/mmterpstra/pipeline-util

=head1 AUTHOR

m.m.terpstra

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2018 by m.m.terpstra

LGPL licence

=cut
