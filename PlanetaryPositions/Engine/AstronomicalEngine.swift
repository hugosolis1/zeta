import Foundation

// MARK: - Astronomical Constants
struct AstroConstants {
    static let J2000: Double = 2451545.0
    static let DEG_TO_RAD: Double = .pi / 180.0
    static let RAD_TO_DEG: Double = 180.0 / .pi
    static let ARCSEC_TO_DEG: Double = 1.0 / 3600.0
}

// MARK: - Planet Data Model
struct PlanetPosition: Identifiable {
    let id = UUID()
    let name: String
    let symbol: String
    let longitude: Double      // ecliptic longitude degrees 0-360
    let latitude: Double       // ecliptic latitude degrees
    let distance: Double       // AU
    let speed: Double          // degrees/day (negative = retrograde)
    let sign: ZodiacSign
    let degreeInSign: Double
    let minuteInSign: Double
    let secondInSign: Double
    var isRetrograde: Bool { speed < 0 }
    
    var formattedPosition: String {
        let d = Int(degreeInSign)
        let m = Int(minuteInSign)
        let s = Int(secondInSign)
        return "\(d)°\(m)'\(s)\" \(sign.rawValue)"
    }
    
    var decimalDegreeInSign: Double { degreeInSign + minuteInSign / 60.0 + secondInSign / 3600.0 }
}

// MARK: - Zodiac
enum ZodiacSign: String, CaseIterable {
    case aries = "♈ Aries"
    case taurus = "♉ Tauro"
    case gemini = "♊ Géminis"
    case cancer = "♋ Cáncer"
    case leo = "♌ Leo"
    case virgo = "♍ Virgo"
    case libra = "♎ Libra"
    case scorpio = "♏ Escorpio"
    case sagittarius = "♐ Sagitario"
    case capricorn = "♑ Capricornio"
    case aquarius = "♒ Acuario"
    case pisces = "♓ Piscis"
    
    var glyph: String {
        switch self {
        case .aries: return "♈"
        case .taurus: return "♉"
        case .gemini: return "♊"
        case .cancer: return "♋"
        case .leo: return "♌"
        case .virgo: return "♍"
        case .libra: return "♎"
        case .scorpio: return "♏"
        case .sagittarius: return "♐"
        case .capricorn: return "♑"
        case .aquarius: return "♒"
        case .pisces: return "♓"
        }
    }
    
    static func from(longitude: Double) -> ZodiacSign {
        let norm = longitude.truncatingRemainder(dividingBy: 360)
        let pos = norm < 0 ? norm + 360 : norm
        let idx = Int(pos / 30.0) % 12
        return ZodiacSign.allCases[idx]
    }
}

// MARK: - Chart Angles
struct ChartAngles {
    let ascendant: Double
    let descendant: Double
    let midheaven: Double
    let imumCoeli: Double
    let northNode: Double
    let southNode: Double
}

// MARK: - Astronomical Engine
class AstronomicalEngine {
    // Constant aliases for use in extensions
    static let DEG_TO_RAD = AstroConstants.DEG_TO_RAD
    static let RAD_TO_DEG = AstroConstants.RAD_TO_DEG
    static let J2000      = AstroConstants.J2000
    
    // Convert calendar date to Julian Day Number
    static func julianDay(date: Date) -> Double {
        let cal = Calendar(identifier: .gregorian)
        let comps = cal.dateComponents([.year, .month, .day, .hour, .minute, .second], from: date)
        
        let Y = Double(comps.year ?? 2000)
        let M = Double(comps.month ?? 1)
        let D = Double(comps.day ?? 1)
        let h = Double(comps.hour ?? 12)
        let m = Double(comps.minute ?? 0)
        let s = Double(comps.second ?? 0)
        
        let dayFraction = D + h/24.0 + m/1440.0 + s/86400.0
        
        var jY = Y
        var jM = M
        if M <= 2 { jY -= 1; jM += 12 }
        
        let A = Int(jY / 100)
        let B = 2 - A + Int(A / 4)
        
        return Double(Int(365.25 * (jY + 4716))) + Double(Int(30.6001 * (jM + 1))) + dayFraction + Double(B) - 1524.5
    }
    
    // Julian centuries from J2000
    static func T(jd: Double) -> Double {
        return (jd - AstroConstants.J2000) / 36525.0
    }
    
    // Normalize angle to 0-360
    static func norm360(_ angle: Double) -> Double {
        var a = angle.truncatingRemainder(dividingBy: 360)
        if a < 0 { a += 360 }
        return a
    }
    
    // Convert decimal degrees to DMS
    static func toDMS(degrees: Double) -> (d: Double, m: Double, s: Double) {
        let total = abs(degrees)
        let d = floor(total)
        let mFrac = (total - d) * 60
        let m = floor(mFrac)
        let s = (mFrac - m) * 60
        return (d, m, s)
    }
    
    // MARK: - Sun (geometric mean, full VSOP87 simplified)
    static func sunLongitude(T: Double) -> Double {
        let L0 = norm360(280.46646 + 36000.76983 * T + 0.0003032 * T * T)
        let M = norm360(357.52911 + 35999.05029 * T - 0.0001537 * T * T) * AstroConstants.DEG_TO_RAD
        let C = (1.914602 - 0.004817 * T - 0.000014 * T * T) * sin(M)
              + (0.019993 - 0.000101 * T) * sin(2 * M)
              + 0.000289 * sin(3 * M)
        
        let sunLon = L0 + C
        // Apparent longitude (aberration)
        let omega = (125.04 - 1934.136 * T) * AstroConstants.DEG_TO_RAD
        let apparent = sunLon - 0.00569 - 0.00478 * sin(omega)
        return norm360(apparent)
    }
    
    // MARK: - Moon (full ELP2000-82 truncated)
    static func moonLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double) {
        let D  = norm360(297.85036 + 445267.111480 * T - 0.0019142 * T*T) * AstroConstants.DEG_TO_RAD
        let M  = norm360(357.52772 + 35999.050340 * T - 0.0001603 * T*T) * AstroConstants.DEG_TO_RAD
        let Mp = norm360(134.96298 + 477198.867398 * T + 0.0086972 * T*T) * AstroConstants.DEG_TO_RAD
        let F  = norm360(93.27191  + 483202.017538 * T - 0.0036825 * T*T) * AstroConstants.DEG_TO_RAD
        
        // Longitude terms (degrees * 10^-4)
        var sumL: Double = 0
        let lTerms: [(Double, Double, Double, Double, Double)] = [
            (0,0,1,0,  6288774),
            (2,0,-1,0, 1274027),
            (2,0,0,0,   658314),
            (0,0,2,0,   213618),
            (0,1,0,0,  -185116),
            (0,0,0,2,  -114332),
            (2,0,-2,0,   58793),
            (2,-1,-1,0,  57066),
            (2,0,1,0,    53322),
            (2,-1,0,0,   45758),
            (0,1,-1,0,  -40923),
            (1,0,0,0,   -34720),
            (0,1,1,0,   -30383),
            (2,0,0,-2,   15327),
            (0,0,1,2,   -12528),
            (0,0,1,-2,   10980),
            (4,0,-1,0,   10675),
            (0,0,3,0,    10034),
            (4,0,-2,0,    8548),
            (2,1,-1,0,   -7888)
        ]
        for t in lTerms {
            let arg = t.0*D + t.1*M + t.2*Mp + t.3*F
            sumL += t.4 * sin(arg)
        }
        
        // Latitude terms
        var sumB: Double = 0
        let bTerms: [(Double, Double, Double, Double, Double)] = [
            (0,0,0,1,  5128122),
            (0,0,1,1,   280602),
            (0,0,1,-1,  277693),
            (2,0,0,-1,  173237),
            (2,0,-1,1,   55413),
            (2,0,-1,-1,  46271),
            (2,0,0,1,   32573),
            (0,0,2,1,   17198),
            (2,0,1,-1,   9266),
            (0,0,2,-1,   8822)
        ]
        for t in bTerms {
            let arg = t.0*D + t.1*M + t.2*Mp + t.3*F
            sumB += t.4 * sin(arg)
        }
        
        // Distance terms
        var sumR: Double = 0
        let rTerms: [(Double, Double, Double, Double, Double)] = [
            (0,0,1,0,  -20905355),
            (2,0,-1,0,  -3699111),
            (2,0,0,0,   -2955968),
            (0,0,2,0,    -569925),
            (0,1,0,0,     48888),
            (0,0,0,2,    -3149),
            (2,0,-2,0,   246158),
            (2,-1,-1,0, -152138),
            (2,0,1,0,   -170733),
            (2,-1,0,0,  -204586)
        ]
        for t in rTerms {
            let arg = t.0*D + t.1*M + t.2*Mp + t.3*F
            sumR += t.4 * cos(arg)
        }
        
        let L1 = norm360(218.3164477 + 481267.88123421 * T)
        let lon = norm360(L1 + sumL / 1000000.0)
        let lat = sumB / 1000000.0
        let dist = (385000.56 + sumR / 1000.0) / 149597870.7 // AU
        
        return (lon, lat, dist)
    }
    
    // MARK: - Outer planets using truncated VSOP87
    // Mercury
    static func mercuryLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(252.250906 + 149474.0722491 * T + 0.0003035 * T*T)
        let a = 0.387098310
        let e = 0.20563175 + 0.000020407 * T
        let M = norm360(174.791086 + 149472.515654 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 4.09)
    }
    
    // Venus
    static func venusLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(181.979801 + 58519.2130302 * T + 0.00031014 * T*T)
        let a = 0.723329820
        let e = 0.00677188 - 0.000047766 * T
        let M = norm360(50.416009 + 58517.803876 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 1.60)
    }
    
    // Mars
    static func marsLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(355.433275 + 19141.6964746 * T + 0.00031052 * T*T)
        let a = 1.523679342
        let e = 0.09340062 + 0.000090483 * T
        let M = norm360(19.373481 + 19140.2993313 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 0.524)
    }
    
    // Jupiter
    static func jupiterLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(34.351519 + 3036.3027748 * T + 0.00022330 * T*T)
        let a = 5.202603191 + 0.0000001913 * T
        let e = 0.04849485 + 0.000163244 * T
        let M = norm360(20.9240022 + 3034.9056606 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 0.0831)
    }
    
    // Saturn
    static func saturnLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(50.077444 + 1223.5110686 * T + 0.00051908 * T*T)
        let a = 9.554909596 - 0.0000021389 * T
        let e = 0.05550862 - 0.000346818 * T
        let M = norm360(317.020831 + 1222.1134582 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 0.0335)
    }
    
    // Uranus
    static func uranusLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(314.055005 + 429.8640561 * T + 0.00030390 * T*T)
        let a = 19.218446062 - 0.0000000372 * T
        let e = 0.04629590 - 0.000027337 * T
        let M = norm360(142.5905632 + 428.4603567 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 0.01176)
    }
    
    // Neptune
    static func neptuneLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let L = norm360(304.348665 + 219.8833092 * T + 0.00030882 * T*T)
        let a = 30.110386869 - 0.0000001663 * T
        let e = 0.00898809 + 0.000006408 * T
        let M = norm360(256.2281046 + 218.4862002 * T) * AstroConstants.DEG_TO_RAD
        let v = trueAnomaly(M: M, e: e)
        let r = a * (1 - e*e) / (1 + e * cos(v))
        let lon = norm360(L + (v - M) * AstroConstants.RAD_TO_DEG)
        return (lon, 0, r, 0.006)
    }
    
    // Pluto (simplified)
    static func plutoLongitude(T: Double) -> (lon: Double, lat: Double, dist: Double, speed: Double) {
        let lon = norm360(238.9508 + 144.9600 * T)
        return (lon, -17.0, 39.543, 0.004)
    }
    
    // True anomaly from mean anomaly (Kepler equation solver)
    static func trueAnomaly(M: Double, e: Double) -> Double {
        var E = M
        for _ in 0..<50 {
            let dE = (M - E + e * sin(E)) / (1 - e * cos(E))
            E += dE
            if abs(dE) < 1e-10 { break }
        }
        let sinV = sqrt(1 - e*e) * sin(E) / (1 - e * cos(E))
        let cosV = (cos(E) - e) / (1 - e * cos(E))
        return atan2(sinV, cosV)
    }
    
    // MARK: - Geocentric conversion
    static func toGeocentric(helioLon: Double, helioLat: Double, helioDist: Double,
                              sunLon: Double, sunDist: Double) -> (lon: Double, lat: Double, dist: Double) {
        let lon1 = helioLon * AstroConstants.DEG_TO_RAD
        let lat1 = helioLat * AstroConstants.DEG_TO_RAD
        let sLon = sunLon * AstroConstants.DEG_TO_RAD
        
        let x = helioDist * cos(lat1) * cos(lon1) - sunDist * cos(sLon)
        let y = helioDist * cos(lat1) * sin(lon1) - sunDist * sin(sLon)
        let z = helioDist * sin(lat1)
        
        let dist = sqrt(x*x + y*y + z*z)
        let lon = norm360(atan2(y, x) * AstroConstants.RAD_TO_DEG)
        let lat = atan2(z, sqrt(x*x + y*y)) * AstroConstants.RAD_TO_DEG
        
        return (lon, lat, dist)
    }
    
    // MARK: - Lunar Nodes
    static func lunarNodes(T: Double) -> (north: Double, south: Double) {
        let omega = norm360(125.0445479 - 1934.1362608 * T + 0.0020754 * T*T + T*T*T/467441.0)
        return (omega, norm360(omega + 180))
    }
    
    // MARK: - Obliquity of Ecliptic
    static func obliquity(T: Double) -> Double {
        return 23.4392911 - 0.0130042 * T - 0.00000164 * T*T + 0.000000504 * T*T*T
    }
    
    // MARK: - GMST (Greenwich Mean Sidereal Time)
    static func gmst(jd: Double) -> Double {
        let T = (jd - AstroConstants.J2000) / 36525.0
        let gmst0 = 280.46061837 + 360.98564736629 * (jd - AstroConstants.J2000)
                  + 0.000387933 * T*T - T*T*T / 38710000.0
        return norm360(gmst0)
    }
    
    // MARK: - Local Sidereal Time
    static func lst(jd: Double, longitude: Double) -> Double {
        return norm360(gmst(jd: jd) + longitude)
    }
    
    // MARK: - Ascendant calculation
    static func ascendant(obliquity: Double, lst: Double, latitude: Double) -> Double {
        let E = obliquity * AstroConstants.DEG_TO_RAD
        let T = lst * AstroConstants.DEG_TO_RAD
        let L = latitude * AstroConstants.DEG_TO_RAD
        
        let y = -cos(T)
        let x = sin(E) * tan(L) + cos(E) * sin(T)
        var asc = atan2(y, x) * AstroConstants.RAD_TO_DEG
        if asc < 0 { asc += 360 }
        return norm360(asc)
    }
    
    // MARK: - Midheaven
    static func midheaven(obliquity: Double, lst: Double) -> Double {
        let E = obliquity * AstroConstants.DEG_TO_RAD
        let T = lst * AstroConstants.DEG_TO_RAD
        let mc = atan2(sin(T), cos(T) * cos(E) - tan(0) * sin(E)) * AstroConstants.RAD_TO_DEG
        return norm360(mc)
    }
    
    // MARK: - Main compute function
    // Uses high-precision algorithms from SwissEphExtension.
    // For geocentric: heliocentric coordinates are converted to geocentric via 3D vector subtraction.
    // For heliocentric: raw heliocentric ecliptic coordinates are returned.
    // The Sun's heliocentric longitude is always the Earth's heliocentric longitude + 180°.
    static func computePositions(date: Date, isGeocentric: Bool,
                                  latitude: Double = 19.4326,
                                  longitude: Double = -99.1332) -> ([PlanetPosition], ChartAngles) {
        let jd = julianDay(date: date)
        let t = T(jd: jd)
        
        // --- Earth/Sun ---
        // Geocentric Sun longitude (apparent, for display and geocentric conversion)
        let geoSunLon = sunLongitudeHighPrecision(T: t)
        // Earth's heliocentric distance (for geocentric conversion)
        let Me_rad = norm360(357.52911 + 35999.05029 * t) * AstroConstants.DEG_TO_RAD
        let earthDist = 1.000001018 * (1 - 0.016708634 * 0.016708634) / (1 + 0.016708634 * cos(Me_rad))
        // Earth's heliocentric longitude = geocentric Sun longitude + 180°
        let earthHelioLon = norm360(geoSunLon + 180.0)
        
        let obliq = obliquity(T: t)
        let lstDeg = lst(jd: jd, longitude: longitude)
        let asc = ascendant(obliquity: obliq, lst: lstDeg, latitude: latitude)
        let mc = midheaven(obliquity: obliq, lst: lstDeg)
        let nodes = lunarNodes(T: t)
        
        // Helper: make PlanetPosition from heliocentric coords, applying geo conversion if needed.
        // Numeric speed is estimated by computing position 1 day later.
        func makePlanet(name: String, symbol: String,
                        hLon: Double, hLat: Double, hDist: Double,
                        hLon2: Double, hLat2: Double, hDist2: Double) -> PlanetPosition {
            let lon: Double
            let lat: Double
            let dist: Double
            let speed: Double
            
            if isGeocentric {
                // Current geocentric position
                let geo = toGeocentric(helioLon: hLon, helioLat: hLat, helioDist: hDist,
                                       sunLon: earthHelioLon, sunDist: earthDist)
                // Position 1 day later for speed (earth also moves ~0.9856°/day)
                let earthHelioLon2 = norm360(earthHelioLon + 0.9856)
                let geo2 = toGeocentric(helioLon: hLon2, helioLat: hLat2, helioDist: hDist2,
                                        sunLon: earthHelioLon2, sunDist: earthDist)
                lon = geo.lon
                lat = geo.lat
                dist = geo.dist
                // Angular speed: signed difference for correct retrograde detection
                var diff = geo2.lon - geo.lon
                if diff > 180 { diff -= 360 }
                if diff < -180 { diff += 360 }
                speed = diff
            } else {
                lon = hLon
                lat = hLat
                dist = hDist
                var diff = hLon2 - hLon
                if diff > 180 { diff -= 360 }
                if diff < -180 { diff += 360 }
                speed = diff
            }
            
            let sign = ZodiacSign.from(longitude: lon)
            let degInSign = lon.truncatingRemainder(dividingBy: 30)
            let dms = toDMS(degrees: degInSign < 0 ? degInSign + 30 : degInSign)
            return PlanetPosition(name: name, symbol: symbol,
                                  longitude: lon, latitude: lat, distance: dist, speed: speed,
                                  sign: sign, degreeInSign: dms.d, minuteInSign: dms.m, secondInSign: dms.s)
        }
        
        // High-precision heliocentric positions for T and T+1day
        let t2 = T(jd: jd + 1.0)
        
        // Sun/Earth: in heliocentric mode show Earth's position (= Sun at 0, Earth orbits Sun)
        // In astrology heliocentric, "Sun" is placed at 0° Aries (centre), so we show Earth instead.
        let sunLon: Double
        let sunLat: Double = 0.0
        let sunDist: Double
        let sunSpeedDisplay: Double
        
        if isGeocentric {
            // Geocentric: Sun appears at geoSunLon
            let geoSunLon2 = sunLongitudeHighPrecision(T: t2)
            var diff = geoSunLon2 - geoSunLon; if diff > 180 { diff -= 360 }; if diff < -180 { diff += 360 }
            sunLon = geoSunLon
            sunDist = earthDist
            sunSpeedDisplay = diff
        } else {
            // Heliocentric: Earth's longitude (= geocentric Sun + 180°)
            sunLon = earthHelioLon
            sunDist = earthDist
            let earthHelioLon2 = norm360(sunLongitudeHighPrecision(T: t2) + 180.0)
            var diff = earthHelioLon2 - earthHelioLon; if diff > 180 { diff -= 360 }; if diff < -180 { diff += 360 }
            sunSpeedDisplay = diff
        }
        let sunSign = ZodiacSign.from(longitude: sunLon)
        let sunDms = toDMS(degrees: sunLon.truncatingRemainder(dividingBy: 30))
        let sunLabel = isGeocentric ? "Sol" : "Tierra"
        let sunSymbol = isGeocentric ? "☉" : "⊕"
        let sunPos = PlanetPosition(name: sunLabel, symbol: sunSymbol,
                                    longitude: sunLon, latitude: sunLat, distance: sunDist, speed: sunSpeedDisplay,
                                    sign: sunSign, degreeInSign: sunDms.d, minuteInSign: sunDms.m, secondInSign: sunDms.s)
        
        // Moon (always geocentric; not shown in heliocentric mode)
        let moonData  = moonLongitudeHighPrecision(T: t)
        let moonData2 = moonLongitudeHighPrecision(T: t2)
        var moonDiff = moonData2.lon - moonData.lon; if moonDiff > 180 { moonDiff -= 360 }; if moonDiff < -180 { moonDiff += 360 }
        let moonSign = ZodiacSign.from(longitude: moonData.lon)
        let moonDms  = toDMS(degrees: moonData.lon.truncatingRemainder(dividingBy: 30))
        let moonPos  = PlanetPosition(name: "Luna", symbol: "☽",
                                      longitude: moonData.lon, latitude: moonData.lat,
                                      distance: moonData.dist * 149597870.7,
                                      speed: moonDiff,
                                      sign: moonSign,
                                      degreeInSign: moonDms.d, minuteInSign: moonDms.m, secondInSign: moonDms.s)
        
        // High-precision heliocentric planet positions
        let merc  = mercuryHighPrecision(T: t);  let merc2  = mercuryHighPrecision(T: t2)
        let ven   = venusHighPrecision(T: t);    let ven2   = venusHighPrecision(T: t2)
        let mars  = marsHighPrecision(T: t);     let mars2  = marsHighPrecision(T: t2)
        let jup   = jupiterHighPrecision(T: t);  let jup2   = jupiterHighPrecision(T: t2)
        let sat   = saturnHighPrecision(T: t);   let sat2   = saturnHighPrecision(T: t2)
        let ura   = uranusHighPrecision(T: t);   let ura2   = uranusHighPrecision(T: t2)
        let nep   = neptuneHighPrecision(T: t);  let nep2   = neptuneHighPrecision(T: t2)
        let plu   = plutoHighPrecision(T: t);    let plu2   = plutoHighPrecision(T: t2)
        
        var planets: [PlanetPosition] = [sunPos]
        
        // Moon only in geocentric mode
        if isGeocentric {
            planets.append(moonPos)
        }
        
        planets += [
            makePlanet(name: "Mercurio", symbol: "☿",
                       hLon: merc.lon, hLat: merc.lat, hDist: merc.dist,
                       hLon2: merc2.lon, hLat2: merc2.lat, hDist2: merc2.dist),
            makePlanet(name: "Venus",    symbol: "♀",
                       hLon: ven.lon,  hLat: ven.lat,  hDist: ven.dist,
                       hLon2: ven2.lon,  hLat2: ven2.lat,  hDist2: ven2.dist),
            makePlanet(name: "Marte",    symbol: "♂",
                       hLon: mars.lon, hLat: mars.lat, hDist: mars.dist,
                       hLon2: mars2.lon, hLat2: mars2.lat, hDist2: mars2.dist),
            makePlanet(name: "Júpiter",  symbol: "♃",
                       hLon: jup.lon,  hLat: jup.lat,  hDist: jup.dist,
                       hLon2: jup2.lon,  hLat2: jup2.lat,  hDist2: jup2.dist),
            makePlanet(name: "Saturno",  symbol: "♄",
                       hLon: sat.lon,  hLat: sat.lat,  hDist: sat.dist,
                       hLon2: sat2.lon,  hLat2: sat2.lat,  hDist2: sat2.dist),
            makePlanet(name: "Urano",    symbol: "♅",
                       hLon: ura.lon,  hLat: ura.lat,  hDist: ura.dist,
                       hLon2: ura2.lon,  hLat2: ura2.lat,  hDist2: ura2.dist),
            makePlanet(name: "Neptuno",  symbol: "♆",
                       hLon: nep.lon,  hLat: nep.lat,  hDist: nep.dist,
                       hLon2: nep2.lon,  hLat2: nep2.lat,  hDist2: nep2.dist),
            makePlanet(name: "Plutón",   symbol: "♇",
                       hLon: plu.lon,  hLat: plu.lat,  hDist: plu.dist,
                       hLon2: plu2.lon,  hLat2: plu2.lat,  hDist2: plu2.dist),
        ]
        
        // Lunar nodes only in geocentric mode
        if isGeocentric {
            let nnDms = toDMS(degrees: nodes.north.truncatingRemainder(dividingBy: 30))
            let nnPos = PlanetPosition(name: "Nodo Norte", symbol: "☊",
                                       longitude: nodes.north, latitude: 0, distance: 0, speed: -0.053,
                                       sign: ZodiacSign.from(longitude: nodes.north),
                                       degreeInSign: nnDms.d, minuteInSign: nnDms.m, secondInSign: nnDms.s)
            let snDms = toDMS(degrees: nodes.south.truncatingRemainder(dividingBy: 30))
            let snPos = PlanetPosition(name: "Nodo Sur", symbol: "☋",
                                       longitude: nodes.south, latitude: 0, distance: 0, speed: -0.053,
                                       sign: ZodiacSign.from(longitude: nodes.south),
                                       degreeInSign: snDms.d, minuteInSign: snDms.m, secondInSign: snDms.s)
            planets.append(nnPos)
            planets.append(snPos)
        }
        
        let angles = ChartAngles(
            ascendant: asc,
            descendant: norm360(asc + 180),
            midheaven: mc,
            imumCoeli: norm360(mc + 180),
            northNode: nodes.north,
            southNode: nodes.south
        )
        
        return (planets, angles)
    }
    
    // MARK: - Placidus House Cusps (Exact Calculation)
    static func placidusHouses(obliquity: Double, lst: Double, latitude: Double) -> [Double] {
        let eps = obliquity * AstroConstants.DEG_TO_RAD
        let phi = latitude * AstroConstants.DEG_TO_RAD
        let ramc = lst * AstroConstants.DEG_TO_RAD
        
        // MC is the 10th house cusp
        let mc = atan2(sin(ramc), cos(ramc) * cos(eps)) * AstroConstants.RAD_TO_DEG
        let mcNorm = norm360(mc)
        
        // ASC is the 1st house cusp
        let asc = ascendant(obliquity: obliquity, lst: lst, latitude: latitude)
        
        var cusps: [Double] = Array(repeating: 0, count: 12)
        cusps[0] = asc   // House 1
        cusps[9] = mcNorm // House 10
        
        // Calculate intermediate house cusps using Placidus method
        // Houses 11, 12 (above horizon)
        for i in 1...2 {
            let factor = Double(i) / 3.0
            var ra = ramc + factor * .pi / 2
            var decl: Double
            var iteration = 0
            
            // Iterative solution for Placidus
            while iteration < 20 {
                let x = cos(ra) * cos(eps)
                decl = asin(sin(ra) * sin(eps))
                
                let y = atan2(tan(phi), cos(decl) * cos(x != 0 ? acos(-x / cos(decl)) : 0))
                let ha = factor * .pi / 2
                
                ra = ramc + ha - y
                iteration += 1
            }
            
            let cusp = atan2(sin(ra), cos(ra) * cos(eps)) * AstroConstants.RAD_TO_DEG
            cusps[9 + i] = norm360(cusp) // Houses 11, 12
        }
        
        // Houses 2, 3 (below horizon)
        for i in 1...2 {
            let factor = Double(i) / 3.0
            var ra = ramc + .pi + factor * .pi / 2
            var iteration = 0
            
            while iteration < 20 {
                let x = cos(ra) * cos(eps)
                let decl = asin(sin(ra) * sin(eps))
                let y = atan2(tan(phi), cos(decl) * cos(x != 0 ? acos(-x / cos(decl)) : 0))
                let ha = factor * .pi / 2
                
                ra = ramc + .pi + ha - y
                iteration += 1
            }
            
            let cusp = atan2(sin(ra), cos(ra) * cos(eps)) * AstroConstants.RAD_TO_DEG
            cusps[i] = norm360(cusp) // Houses 2, 3
        }
        
        // Calculate remaining houses (oppositions)
        cusps[6] = norm360(cusps[0] + 180)  // House 7 opposite House 1
        cusps[7] = norm360(cusps[1] + 180)  // House 8 opposite House 2
        cusps[8] = norm360(cusps[2] + 180)  // House 9 opposite House 3
        cusps[3] = norm360(cusps[9] + 180)  // House 4 opposite House 10
        cusps[4] = norm360(cusps[10] + 180) // House 5 opposite House 11
        cusps[5] = norm360(cusps[11] + 180) // House 6 opposite House 12
        
        return cusps.map { norm360($0) }
    }
    
    // Simplified but accurate house cusps
    static func calculateHouses(date: Date, latitude: Double, longitude: Double) -> [Double] {
        let jd = julianDay(date: date)
        let t = T(jd: jd)
        let obliq = obliquity(T: t)
        let lstDeg = lst(jd: jd, longitude: longitude)
        
        return placidusHouses(obliquity: obliq, lst: lstDeg, latitude: latitude)
    }
    
    // MARK: - Planetary Aspects
    enum AspectType: String {
        case conjunction = "Conjunción"
        case opposition = "Oposición"
        case trine = "Trino"
        case square = "Cuadratura"
        case sextile = "Sextil"
        case quincunx = "Quincuncio"
        case semisextile = "Semisextil"
        case semisquare = "Semicuadratura"
        
        var symbol: String {
            switch self {
            case .conjunction: return "☌"
            case .opposition: return "☍"
            case .trine: return "△"
            case .square: return "□"
            case .sextile: return "✶"
            case .quincunx: return "⚹"
            case .semisextile: return "⚶"
            case .semisquare: return "∠"
            }
        }
        
        var angle: Double {
            switch self {
            case .conjunction: return 0
            case .semisextile: return 30
            case .semisquare: return 45
            case .sextile: return 60
            case .square: return 90
            case .trine: return 120
            case .quincunx: return 150
            case .opposition: return 180
            }
        }
        
        var orb: Double {
            switch self {
            case .conjunction: return 10
            case .opposition: return 8
            case .trine: return 8
            case .square: return 8
            case .sextile: return 6
            case .quincunx: return 3
            case .semisextile: return 3
            case .semisquare: return 2
            }
        }
        
        var color: String {
            switch self {
            case .conjunction: return "green"
            case .opposition: return "red"
            case .trine: return "blue"
            case .square: return "orange"
            case .sextile: return "cyan"
            case .quincunx: return "purple"
            case .semisextile: return "gray"
            case .semisquare: return "brown"
            }
        }
    }
    
    struct Aspect: Identifiable {
        let id = UUID()
        let planet1: String
        let planet1Symbol: String
        let planet2: String
        let planet2Symbol: String
        let aspect: AspectType
        let angle: Double
        let orb: Double
        let applying: Bool
    }
    
    // Calculate angular difference between two longitudes
    static func angularDifference(_ lon1: Double, _ lon2: Double) -> Double {
        var diff = abs(lon1 - lon2)
        if diff > 180 { diff = 360 - diff }
        return diff
    }
    
    // Check if aspect is applying (planets moving closer) or separating
    static func isApplying(lon1: Double, speed1: Double, lon2: Double, speed2: Double, exactAngle: Double) -> Bool {
        let diff = lon2 - lon1
        var normDiff = diff
        while normDiff > 180 { normDiff -= 360 }
        while normDiff < -180 { normDiff += 360 }
        
        let relativeSpeed = speed2 - speed1
        return (normDiff > 0 && relativeSpeed < 0) || (normDiff < 0 && relativeSpeed > 0)
    }
    
    // Calculate all aspects between planets
    static func calculateAspects(planets: [PlanetPosition]) -> [Aspect] {
        var aspects: [Aspect] = []
        
        for i in 0..<planets.count {
            for j in (i+1)..<planets.count {
                let p1 = planets[i]
                let p2 = planets[j]
                
                let diff = angularDifference(p1.longitude, p2.longitude)
                
                for aspectType in [AspectType.conjunction, AspectType.sextile, AspectType.square, AspectType.trine, AspectType.opposition, AspectType.quincunx] {
                    let exactAngle = aspectType.angle
                    let orb = aspectType.orb
                    let deviation = abs(diff - exactAngle)
                    
                    if deviation <= orb {
                        let applying = isApplying(
                            lon1: p1.longitude, speed1: p1.speed,
                            lon2: p2.longitude, speed2: p2.speed,
                            exactAngle: exactAngle
                        )
                        
                        aspects.append(Aspect(
                            planet1: p1.name,
                            planet1Symbol: p1.symbol,
                            planet2: p2.name,
                            planet2Symbol: p2.symbol,
                            aspect: aspectType,
                            angle: exactAngle,
                            orb: deviation,
                            applying: applying
                        ))
                        break // Only one aspect per pair
                    }
                }
            }
        }
        
        // Sort by orb (tightest aspects first)
        return aspects.sorted { $0.orb < $1.orb }
    }
}
