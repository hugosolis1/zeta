import SwiftUI
import Combine

class AstroViewModel: ObservableObject {
    @Published var selectedDate: Date = Date()
    @Published var isGeocentric: Bool = true
    @Published var latitude: Double = 19.4326   // default CDMX
    @Published var longitude: Double = -99.1332
    @Published var planets: [PlanetPosition] = []
    @Published var angles: ChartAngles = ChartAngles(ascendant: 0, descendant: 180, midheaven: 90, imumCoeli: 270, northNode: 0, southNode: 180)
    @Published var houseCusps: [Double] = Array(repeating: 0, count: 12)
    @Published var aspects: [AstronomicalEngine.Aspect] = []
    @Published var isLoading: Bool = false
    
    init() {
        compute()
    }
    
    func compute() {
        isLoading = true
        DispatchQueue.global(qos: .userInitiated).async {
            let result = AstronomicalEngine.computePositions(
                date: self.selectedDate,
                isGeocentric: self.isGeocentric,
                latitude: self.latitude,
                longitude: self.longitude
            )
            
            let houses = AstronomicalEngine.calculateHouses(
                date: self.selectedDate,
                latitude: self.latitude,
                longitude: self.longitude
            )
            
            let aspects = AstronomicalEngine.calculateAspects(planets: result.0)
            
            DispatchQueue.main.async {
                self.planets = result.0
                self.angles = result.1
                self.houseCusps = houses
                self.aspects = aspects
                self.isLoading = false
            }
        }
    }
    
    func anglePosition(name: String, symbol: String, lon: Double) -> PlanetPosition {
        let sign = ZodiacSign.from(longitude: lon)
        let deg = lon.truncatingRemainder(dividingBy: 30)
        let dms = AstronomicalEngine.toDMS(degrees: deg)
        return PlanetPosition(name: name, symbol: symbol,
                              longitude: lon, latitude: 0, distance: 0, speed: 0,
                              sign: sign, degreeInSign: dms.d, minuteInSign: dms.m, secondInSign: dms.s)
    }
    
    var ascendantPosition: PlanetPosition {
        anglePosition(name: "Ascendente", symbol: "AC", lon: angles.ascendant)
    }
    var descendantPosition: PlanetPosition {
        anglePosition(name: "Descendente", symbol: "DC", lon: angles.descendant)
    }
    var midheavenPosition: PlanetPosition {
        anglePosition(name: "Medio Cielo", symbol: "MC", lon: angles.midheaven)
    }
    var imumCoeliPosition: PlanetPosition {
        anglePosition(name: "Fondo Cielo", symbol: "IC", lon: angles.imumCoeli)
    }
    var northNodePosition: PlanetPosition {
        anglePosition(name: "Nodo Norte", symbol: "☊", lon: angles.northNode)
    }
    var southNodePosition: PlanetPosition {
        anglePosition(name: "Nodo Sur", symbol: "☋", lon: angles.southNode)
    }
    
    var allPositions: [PlanetPosition] {
        if isGeocentric {
            return planets + [ascendantPosition, descendantPosition, midheavenPosition, imumCoeliPosition, northNodePosition, southNodePosition]
        } else {
            return planets
        }
    }
    
    // Get house for a given longitude
    func getHouse(for longitude: Double) -> Int {
        for i in 0..<12 {
            let cusp1 = houseCusps[i]
            let cusp2 = houseCusps[(i + 1) % 12]
            
            if cusp1 < cusp2 {
                if longitude >= cusp1 && longitude < cusp2 {
                    return i + 1
                }
            } else {
                if longitude >= cusp1 || longitude < cusp2 {
                    return i + 1
                }
            }
        }
        return 1
    }
    
    // Format longitude in 360 degree format
    static func format360(_ longitude: Double) -> String {
        return String(format: "%.2f°", longitude)
    }
}
