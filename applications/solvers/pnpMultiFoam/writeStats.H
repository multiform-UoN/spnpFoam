forAll(fluidRegions,i)
{
  const scalar vol(Foam::gSum(fluidRegions[i].V()));
  Info<< "Fluid region " << fluidRegions[i].name()
      << " Volume " << vol << endl;
  forAll(speciesSolid[i],sp)
  {
    Info << "    Average " << speciesFluid[i][sp].name() << " = "
    << Foam::gSum(cPtrLFluid[i][sp]().field()*fluidRegions[i].V()) / vol
    << endl;
  }
}



forAll(solidRegions,i)
{
  const scalar vol(Foam::gSum(solidRegions[i].V()));
  Info<< "Solid region " << solidRegions[i].name()
      << " Volume " << vol << endl;
  forAll(speciesSolid[i],sp)
  {
    Info << "    Average " << speciesSolid[i][sp].name() << " = "
    << Foam::gSum(cPtrLSolid[i][sp]().field()*solidRegions[i].V()) / vol
    << endl;
  }
}
