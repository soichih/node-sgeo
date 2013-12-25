
//eventuall, I'd like to have similar function listed on google maps api
//https://developers.google.com/maps/documentation/javascript/reference?hl=en

function toRad(angle) {
  return angle * (Math.PI / 180);
}
function toDeg(angle) {
  return angle * 180 / Math.PI;
}

//http://www.movable-type.co.uk/scripts/latlong.html
//takes pair of latlng, and returns distance in kilometers
exports.distanceHavershine = function(p1, p2) {
    var lat1 = p1.lat;
    var lon1 = p1.lng;
    var lat2 = p2.lat;
    var lon2 = p2.lng;

    var R = 6371; // km
    var dLat = toRad(lat2-lat1);
    var dLon = toRad(lon2-lon1);
    var lat1 = toRad(lat1);
    var lat2 = toRad(lat2);

    var a = Math.sin(dLat/2) * Math.sin(dLat/2) +
            Math.sin(dLon/2) * Math.sin(dLon/2) * Math.cos(lat1) * Math.cos(lat2); 
    var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a)); 
    var d = R * c;
    return d;
}

exports.distance = exports.distanceHavershine;

//http://www.movable-type.co.uk/scripts/latlong-vincenty.html
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2012             */
/*                                                                                                */
/* from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the */
/*       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    */
/*       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf                                             */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

/**
 * Calculates geodetic distance between two points specified by latitude/longitude using 
 * Vincenty inverse formula for ellipsoids
 *
 * @param   {Number} lat1, lon1: first point in decimal degrees
 * @param   {Number} lat2, lon2: second point in decimal degrees
 * @returns (Number} distance in metres between points
 */
exports.distanceVincenty = function(p1, p2) {
    var lat1 = p1.lat;
    var lon1 = p1.lng;
    var lat2 = p2.lat;
    var lon2 = p2.lng;

    var a = 6378137, b = 6356752.314245,  f = 1/298.257223563;  // WGS-84 ellipsoid params
    var L = toRad(lon2-lon1);
    var U1 = Math.atan((1-f) * Math.tan(toRad(lat1)));
    var U2 = Math.atan((1-f) * Math.tan(toRad(lat2)));
    var sinU1 = Math.sin(U1), cosU1 = Math.cos(U1);
    var sinU2 = Math.sin(U2), cosU2 = Math.cos(U2);

    var lambda = L, lambdaP, iterLimit = 100;
    do {
        var sinLambda = Math.sin(lambda), cosLambda = Math.cos(lambda);
        var sinSigma = Math.sqrt((cosU2*sinLambda) * (cosU2*sinLambda) + 
          (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda));
        if (sinSigma==0) return 0;  // co-incident points
        var cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        var sigma = Math.atan2(sinSigma, cosSigma);
        var sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        var cosSqAlpha = 1 - sinAlpha*sinAlpha;
        var cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
        if (isNaN(cos2SigmaM)) cos2SigmaM = 0;  // equatorial line: cosSqAlpha=0 (§6)
        var C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1-C) * f * sinAlpha *
          (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
    } while (Math.abs(lambda-lambdaP) > 1e-12 && --iterLimit>0);

    if (iterLimit==0) return NaN  // formula failed to converge

    var uSq = cosSqAlpha * (a*a - b*b) / (b*b);
    var A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
    var B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
    var deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
    B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
    var s = b*A*(sigma-deltaSigma);

    //s = s.toFixed(3); // round to 1mm precision
    return s/1000; //return in km

    /*
    // note: to return initial/final bearings in addition to distance, use something like:
    var fwdAz = Math.atan2(cosU2*sinLambda,  cosU1*sinU2-sinU1*cosU2*cosLambda);
    var revAz = Math.atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
    return { distance: s, initialBearing: fwdAz.toDeg(), finalBearing: revAz.toDeg() };
    */
}

/*
//http://stackoverflow.com/questions/7676845/how-to-find-a-geographic-point-between-two-other-points
lat1/lon1 -- latitude and longitude of point 1 (in degree)
lat2/lon2 -- latitude and longitude of point 2 (in degree)
t - interpolated point witin the path (through 0 - 1.0) defaulted to 0.5
*/
exports.interpolate = function(p1, p2, t) {
    if(t === undefined) {
        t = 0.5;
    }
    var lat1 = p1.lat;
    var lon1 = p1.lng;
    var lat2 = p2.lat;
    var lon2 = p2.lng;

    // Convert microdegrees to radians
    var alatRad=toRad(lat1);
    var alonRad=toRad(lon1);
    var blatRad=toRad(lat2);
    var blonRad=toRad(lon2);

    console.log(alatRad);
    console.log(alonRad);
    console.log(blatRad);
    console.log(blonRad);

    // Calculate distance in longitude
    var dlon = blonRad - alonRad;

    // Calculate common variables
    var alatRadSin = Math.sin(alatRad);
    var blatRadSin = Math.sin(blatRad);
    var alatRadCos = Math.cos(alatRad);
    var blatRadCos = Math.cos(blatRad);
    var dlonCos = Math.cos(dlon);

    // Find distance from A to B
    var distance=Math.acos(alatRadSin*blatRadSin + alatRadCos*blatRadCos * dlonCos);

    // Find bearing from A to B
    var bearing=Math.atan2( Math.sin(dlon) * blatRadCos, alatRadCos*blatRadSin - alatRadSin*blatRadCos*dlonCos);

    // Find new point
    var angularDistance=distance*t;
    var angDistSin=Math.sin(angularDistance);
    var angDistCos=Math.cos(angularDistance);
    var xlatRad = Math.asin( alatRadSin*angDistCos + alatRadCos*angDistSin*Math.cos(bearing) );
    var xlonRad = alonRad + Math.atan2( Math.sin(bearing)*angDistSin*alatRadCos, angDistCos-alatRadSin*Math.sin(xlatRad));

    // Convert radians to deg
    var xlat=toDeg(xlatRad);
    var xlon=toDeg(xlonRad);

    //normalize
    if(xlat>90)xlat=90;
    if(xlat<-90)xlat=-90;
    while(xlon>180)xlon-=360;
    while(xlon<=-180)xlon+=360;
    return {lat:xlat,lng:xlon};
}
