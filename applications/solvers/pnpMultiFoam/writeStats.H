forAll(fluidRegions,i)
{
  forAll(speciesSolid[i],sp)
  {
    Info << "    Concentration of " << speciesFluid[i][sp].name()
    << " in " << fluidRegions[i].name() << " = "
    << speciesTotConcFluid[i][sp] << endl;
  }
}



forAll(solidRegions,i)
{
  forAll(speciesSolid[i],sp)
  {
    Info << "    Concentration " << speciesSolid[i][sp].name()
    << " in " << fluidRegions[i].name() << " = "
    << speciesTotConcSolid[i][sp] << endl;
  }
}
