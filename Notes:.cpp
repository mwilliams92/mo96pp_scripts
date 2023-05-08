Notes:

The outout tree would have a branch for all orruba events
There would need to be a class called TOrruba that takes the orruba branch address as an arguement 
There needs to be a function called GetOrrubaHit() and GetOrrubaMultiplicity() within the TOrruba class, can also have GetSX3Multiplicity(), GetQQQ5Multiplicity() etc
Thee result of these functions are assigne as a TOrrubaHit object
So there needs to be a class called TOrrubaHit.
There can also be classes called TSX3Hit, TQQQ5Hit, TBB10Hit, TBB10StripHit TQQQ5RingHit, TQQQ5SectorHit, TSX3FrontHit, TSX3BackHit
These Hit classes will contain functions like GetEnergy(), GetADC(), GetDetector()

The there needs to be a GetOrrubaHit function in the branch?

class TOrruba {
public:
  int detnum;
  double energy;

  // Define a member function that returns the energy value
  double GetEnergy() const {
    return energy;
  }
};