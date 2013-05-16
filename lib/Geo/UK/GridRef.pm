package Geo::UK::GridRef;
use strict;
use warnings;

use Sub::Exporter -setup => {
    exports => [ 'grid_ref_to_lat_lon', 'lat_lon_to_grid_ref' ],
    groups => {
        default => [ 'grid_ref_to_lat_lon', 'lat_lon_to_grid_ref' ],
    },
};

use constant PI => 4 * atan2(1, 1);

sub _rad_to_deg { $_[0] * 180 / PI }
sub _deg_to_rad { $_[0] * PI / 180 }

use constant {
    # Airy 1830 major and minor semi-axes
    A => 6377563.396,
    B => 6356256.910,

    # NatGrid scale factor on central meridian
    F0 => 0.9996012717,

    # NatGrid true origin is 49ºN,2ºW
    LAT0 => _deg_to_rad(49),
    LON0 => _deg_to_rad(-2),

    # northing & easting of true origin, metres
    N0 => -100000,
    E0 => 400000,
};

use constant {
    # eccentricity squared
    E2 => 1 - B**2 / A**2,

    N => (A - B) / (A + B),
};

use constant {
    CA => 1 + N + 1.25 * N**2 + 1.25 * N**3,
    CB => 3 * N + 3 * N**2 + 2.625 * N**3,
    CC => 1.825 * N**2 + 1.825 * N**3,
    CD => 35 / 24 * N**3,

    C => B * F0,
};

sub lat_lon_to_grid_ref {
    my ($lat, $lon) = @_;

    $lat = _deg_to_rad($lat);
    $lon = _deg_to_rad($lon);

    my $cosLat = cos($lat);
    my $sinLat = sin($lat);
    my $tanLat = $sinLat / $cosLat;

    # transverse radius of curvature
    my $nu = A * F0 / sqrt(1 - E2 * $sinLat**2);

    # meridional radius of curvature
    my $rho = A * F0 * (1 - E2) / sqrt((1 - E2 * $sinLat**2)**3);

    my $eta2 = $nu / $rho - 1;

    my $l1 = $lat - LAT0;
    my $l2 = $lat + LAT0;

    my $Ma = CA * $l1;
    my $Mb = CB * sin($l1) * cos($l2);
    my $Mc = CC * sin(2 * $l1) * cos(2 * $l2);
    my $Md = CD * sin(3 * $l1) * cos(3 * $l2);

    # meridional arc
    my $M = C * ($Ma - $Mb + $Mc - $Md);

    my $I    = $M + N0;
    my $II   = ($nu / 2) * $sinLat * $cosLat;
    my $III  = ($nu / 24) * $sinLat * $cosLat**3 * (5 - $tanLat**2 + 9 * $eta2);
    my $IIIA = ($nu / 720) * $sinLat * $cosLat**5 * (61 - 58 * $tanLat**2 + $tanLat**4);
    my $IV   = $nu * $cosLat;
    my $V    = ($nu / 6) * $cosLat**3 * ($nu / $rho - $tanLat**2);
    my $VI   = ($nu / 120) * $cosLat**5 * (5 - 18 * $tanLat**2 + $tanLat**4 + 14 * $eta2 - 58 * $tanLat**2 * $eta2);

    my $dLon = $lon - LON0;

    my $east  = E0 + $IV * $dLon + $V * $dLon**3 + $VI * $dLon**5;
    my $north = $I + $II * $dLon**2 + $III * $dLon**4 + $IIIA * $dLon**6;

    return _numeric_to_standard($east, $north);
}

sub grid_ref_to_lat_lon {
    my ($grid_ref) = @_;

    my ($east, $north) = _standard_to_numeric($grid_ref);

    my $lat = LAT0;
    my $M   = 0;

    do {
        $lat = ($north - N0 - $M) / (A * F0) + $lat;

        my $l1 = $lat - LAT0;
        my $l2 = $lat + LAT0;

        my $Ma = CA * $l1;
        my $Mb = CB * sin($l1) * cos($l2);
        my $Mc = CC * sin(2 * $l1) * cos(2 * $l2);
        my $Md = CD * sin(3 * $l1) * cos(3 * $l2);

        # meridional arc
        $M = C * ($Ma - $Mb + $Mc - $Md);
    } while ($north - N0 - $M >= 0.00001); # ie until < 0.01mm

    my $cosLat = cos($lat);
    my $sinLat = sin($lat);
    my $tanLat = sin($lat) / cos($lat);
    my $secLat = 1 / $cosLat;

    # transverse radius of curvature
    my $nu = A * F0 / sqrt(1 - E2 * $sinLat**2);

    # meridional radius of curvature
    my $rho = A * F0 * (1 - E2) / sqrt((1 - E2 * $sinLat**2)**3);

    my $eta2 = $nu / $rho - 1;

    my $VII  = $tanLat / (2 * $rho * $nu);
    my $VIII = $tanLat / (24 * $rho * $nu**3) * (5 + 3 * $tanLat**2 + $eta2 - 9 * $tanLat**2 * $eta2);
    my $IX   = $tanLat / (720 * $rho * $nu**5) * (61 + 90 * $tanLat**2 + 45 * $tanLat**4);
    my $X    = $secLat / $nu;
    my $XI   = $secLat / (6 * $nu**3) * ($nu / $rho + 2 * $tanLat**2);
    my $XII  = $secLat / (120 * $nu**5) * (5 + 28 * $tanLat**2 + 24 * $tanLat**4);
    my $XIIA = $secLat / (5040 * $nu**7) * (61 + 662 * $tanLat**2 + 1320 * $tanLat**4 + 720 * $tanLat**6);

    my $dE  = ($east - E0);

    $lat = $lat - $VII * $dE**2 + $VIII * $dE**4 - $IX * $dE**6;
    my $lon = LON0 + $X * $dE - $XI * $dE**3 + $XII * $dE**5 - $XIIA * $dE**7;

    return (_rad_to_deg($lat), _rad_to_deg($lon));
}

sub _standard_to_numeric {
    my ($standard) = @_;

    $standard =~ s/\s+//g;

    # get numeric values of letter references, mapping A->0, B->1, C->2, etc:
    my $l1 = ord(uc(substr($standard, 0, 1))) - ord('A');
    my $l2 = ord(uc(substr($standard, 1, 1))) - ord('A');
    # shuffle down letters after 'I' since 'I' is not used in grid:
    $l1-- if $l1 > 7;
    $l2-- if $l2 > 7;

    # convert grid letters into 100km-square indexes from false origin (grid
    # square SV):
    my $e = (($l1 - 2) % 5) * 5 + ($l2 % 5);
    my $n = (19 - int($l1 / 5) * 5) - int($l2 / 5);

    die "invalid grid reference $standard"
        if $e < 0 || $e > 6 || $n < 0 || $n > 12;

    my $len = length($standard) - 2;

    # append numeric part of references to grid index:
    $e .= substr($standard, 2, $len / 2);
    $n .= substr($standard, 2 + $len / 2);

    # normalise to 1m grid, rounding up to centre of grid square:
    if ($len == 0) {
        $e .= '50000';
        $n .= '50000';
    }
    elsif ($len == 2) {
        $e .= '5000';
        $n .= '5000';
    }
    elsif ($len == 4) {
        $e .= '500';
        $n .= '500';
    }
    elsif ($len == 6) {
        $e .= '50';
        $n .= '50';
    }
    elsif ($len == 8) {
        $e .= '5';
        $n .= '5';
    }
    elsif ($len == 10) {
        # 10-digit refs are already 1m
    }
    else {
        die "invalid grid reference $standard";
    }

    return ($e, $n);
}

sub _numeric_to_standard {
    my ($e, $n, $precision) = @_;
    $precision = 10 unless defined $precision;

    # get the 100km-grid indices
    my $e100k = int($e / 100000);
    my $n100k = int($n / 100000);

    die "invalid position ($e, $n)"
        if $e100k < 0 || $e100k > 6 || $n100k < 0 || $n100k > 12;

    # translate those into numeric equivalents of the grid letters
    my $l1 = (19 - $n100k) - (19 - $n100k) % 5 + int(($e100k + 10) / 5);
    my $l2 = (19 - $n100k) * 5 % 25 + $e100k % 5;

    # compensate for skipped 'I' and build grid letter-pairs
    $l1++ if $l1 > 7;
    $l2++ if $l2 > 7;
    my $letPair = chr($l1 + ord('A')) . chr($l2 + ord('A'));

    my $width = $precision / 2;

    # strip 100km-grid indices from easting & northing, and reduce precision
    $e = int(($e % 100000) / 10**(5 - $width));
    $n = int(($n % 100000) / 10**(5 - $width));

    return sprintf("%s %0${width}d %0${width}d", $letPair, $e, $n);
}

1;
