package proyecto
  model solution "Chemical solution as homogenous mixture of the substances"
    extends Chemical.Icons.Solution;
    extends Chemical.Interfaces.PartialSolution(T(start = temperature_start), p(start = BasePressure));
    parameter Boolean useMechanicPorts = false "Are mechanic ports pressent?" annotation(
      Evaluate = true,
      HideResult = true,
      choices(checkbox = true),
      Dialog(group = "Conditional inputs"));
    parameter Modelica.SIunits.Area SurfaceArea = 0.01 "Area for surfacePort to connect MultiBody components" annotation(
      HideResult = true,
      Dialog(enable = useMechanicPorts));
    parameter Boolean isPistonPositionAbsolute = false "Relavite position has zero at initial state without force" annotation(
      HideResult = true,
      Dialog(enable = useMechanicPorts));
    parameter Boolean useThermalPort = false "Is thermal port pressent?" annotation(
      Evaluate = true,
      HideResult = true,
      choices(checkbox = true),
      Dialog(group = "Conditional inputs"));
    parameter Boolean ConstantTemperature = true "Has the solution constant temperature during simulation (if not useThermalPort)?" annotation(
      HideResult = true,
      Dialog(enable = not useThermalPort));
    parameter Modelica.SIunits.Pressure BasePressure = 100000 "Ambient pressure if useMechanicPort, start pressure or absolute pressure if ConstantPressure" annotation(
      HideResult = true);
    parameter Modelica.SIunits.Temperature temperature_start = 298.15 "Initial temperature of the solution" annotation(
      Dialog(group = "Initialization"));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(T = T, Q_flow = heatFromEnvironment) if useThermalPort annotation(
      Placement(transformation(extent = {{-70, -90}, {-50, -70}}), iconTransformation(extent = {{-62, -104}, {-58, -100}})));
  protected
    parameter Modelica.SIunits.Position positionShift(fixed = false) "=0 absolute, otherwise negative";
    Modelica.SIunits.Position top_s, ds;
    Modelica.SIunits.Force f;
  initial equation
    T = temperature_start;
    positionShift = if isPistonPositionAbsolute then 0 else volume / SurfaceArea;
//s=volume/SurfaceArea - positionShift;
  equation
//hydraulic
    ds = volume / SurfaceArea - positionShift;
    workFromEnvironment = -der(f * ds);
//=der( (p-p0) * volume)
    solution.p = BasePressure - f / SurfaceArea;
    if not useMechanicPorts then
      f = 0;
      top_s = ds;
//equivalent for bottom_s==0
    end if;
//thermal
    if not useThermalPort and ConstantTemperature then
//Ideal thermal exchange between environment and solution to reach constant temperature
      der(T) = 0;
    end if;
    if not useThermalPort and not ConstantTemperature then
//Thermally isolated without any thermal exchange with environment
      heatFromEnvironment = 0;
    end if;
    annotation(
      Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 1, extent = {{-100, -100}, {100, 100}}), graphics = {Text(extent = {{-92, -86}, {76, -94}}, lineColor = {0, 0, 255}, textString = "%name", horizontalAlignment = TextAlignment.Left)}),
      Documentation(revisions = "<html>
  <p>2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
  </html>", info = "<html>
  <h4>amountOfSolution = &sum; amountOfSubstances</h4>
  <h4>mass = &sum; massOfSubstances</h4>
  <h4>volume = &sum; volumeOfSubstances</h4>
  <h4>freeGibbsEnergy = &sum; freeGibbsEnergiesOfSubstances</h4>
  <p>To calculate the sum of extensive substance's properties is misused the Modelica \"flow\" prefix even there are not real physical flows. </p>
  </html>"));
  end solution;

  model substance "Substance in solution"
    extends Chemical.Icons.Substance;
  
    Modelica.SIunits.Concentration c "Molar concentration";
  
    extends Chemical.Interfaces.PartialSubstanceInSolution;
  
    //If it is selected the amount of solution per one kilogram of solvent then the values of amountOfSubstance will be the same as molality
    //If it is selected the amount of solution in one liter of solution then the values of amountOfSubstance will be the same as molarity
    parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start=1e-8
    "Initial amount of the substance in compartment"   annotation(HideResult=true);
  
  protected
    Modelica.SIunits.AmountOfSubstance amountOfSubstance(start=amountOfSubstance_start);
    Real log10n(stateSelect=StateSelect.prefer, start=log10(amountOfSubstance_start))
    "Decadic logarithm of the amount of the substance in solution";
    constant Real InvLog_10=1/log(10);
  
  initial equation
  
    amountOfSubstance=amountOfSubstance_start;
  
  equation
  
    //The main accumulation equation is "der(amountOfSubstance)=port_a.q"
    // However, the numerical solvers can handle it in form of log10n much better. :-)
    der(log10n)=(InvLog_10)*(port_a.q/amountOfSubstance);
    amountOfSubstance = 10^log10n;
  
    //Molar Concentration
    c = amountOfSubstance/solution.V;
  
    //Mole fraction is an analogy of molar concentration or molality.
    x = amountOfSubstance/solution.n;
  
    //solution flows
    solution.dH = molarEnthalpy*port_a.q + der(molarEnthalpy)*amountOfSubstance;
    solution.i = Modelica.Constants.F * z * port_a.q + Modelica.Constants.F*der(z)*amountOfSubstance;
    solution.dV = molarVolume * port_a.q + der(molarVolume)*amountOfSubstance;
  
    //extensive properties
    solution.nj=amountOfSubstance;
    solution.mj=amountOfSubstance*molarMass;
    solution.Vj=amountOfSubstance*molarVolume;
    solution.Gj=amountOfSubstance*port_a.u;
    solution.Qj=Modelica.Constants.F*amountOfSubstance*z;
    solution.Ij=(1/2) * ( amountOfSubstance * z^2);
    solution.otherPropertiesOfSubstance=amountOfSubstance * otherPropertiesPerSubstance;
  
                                                                                                      annotation (
      Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={0,0,255},
          textString="%name")}),
      Documentation(revisions="<html>
  <p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
  </html>", info="<html>
  <h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
  <p>where n is amount of the substance and x is mole fraction.</p>
  <p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
  <p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
  <p>x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
  <p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
  <p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of x will be the same as molality.</p>
  <p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of x will be the same as molarity.</p>
  <p><br><br>Definition of electro-chemical potential:</p>
  <h4>u = u&deg; + R*T*ln(gamma*x) + z*F*v</h4>
  <h4>u&deg; = DfG = DfH - T * DfS</h4>
  <p>where</p>
  <p>x .. mole fraction of the substance in the solution</p>
  <p>T .. temperature in Kelvins</p>
  <p>v .. relative eletric potential of the solution</p>
  <p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
  <p>R .. gas constant</p>
  <p>F .. Faraday constant</p>
  <p>gamma .. activity coefficient</p>
  <p>u&deg; .. chemical potential of pure substance</p>
  <p>DfG .. free Gibbs energy of formation of the substance</p>
  <p>DfH .. free enthalpy of formation of the substance</p>
  <p>DfS .. free entropy of formation of the substance </p>
  <p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
  </html>"));
  end substance;

  model substanceb "Substance in solution"
    extends Chemical.Icons.Substance;
  
    Modelica.SIunits.Concentration c "Molar concentration";
  
    extends Chemical.Interfaces.PartialSubstanceInSolution;
  
    //If it is selected the amount of solution per one kilogram of solvent then the values of amountOfSubstance will be the same as molality
    //If it is selected the amount of solution in one liter of solution then the values of amountOfSubstance will be the same as molarity
    parameter Modelica.SIunits.AmountOfSubstance amountOfSubstance_start=1e-8
    "Initial amount of the substance in compartment"   annotation(HideResult=true);
  
  protected
    Modelica.SIunits.AmountOfSubstance amountOfSubstance(start=amountOfSubstance_start);
    Real log10n(stateSelect=StateSelect.prefer, start=log10(amountOfSubstance_start))
    "Decadic logarithm of the amount of the substance in solution";
    constant Real InvLog_10=1/log(10);
  
  initial equation
  
    amountOfSubstance=amountOfSubstance_start;
  
  equation
  
    //The main accumulation equation is "der(amountOfSubstance)=port_a.q"
    // However, the numerical solvers can handle it in form of log10n much better. :-)
    der(log10n)=(InvLog_10)*(port_a.q/amountOfSubstance);
    amountOfSubstance = 10^log10n;
  
    //Molar Concentration
    c = amountOfSubstance/solution.V;
  
    //Mole fraction is an analogy of molar concentration or molality.
    x = amountOfSubstance/solution.n;
  
    //solution flows
    solution.dH = molarEnthalpy*port_a.q + der(molarEnthalpy)*amountOfSubstance;
    solution.i = Modelica.Constants.F * z * port_a.q + Modelica.Constants.F*der(z)*amountOfSubstance;
    solution.dV = molarVolume * port_a.q + der(molarVolume)*amountOfSubstance;
  
    //extensive properties
    solution.nj=amountOfSubstance;
    solution.mj=amountOfSubstance*molarMass;
    solution.Vj=amountOfSubstance*molarVolume;
    solution.Gj=amountOfSubstance*port_a.u;
    solution.Qj=Modelica.Constants.F*amountOfSubstance*z;
    solution.Ij=(1/2) * ( amountOfSubstance * z^2);
    solution.otherPropertiesOfSubstance=amountOfSubstance * otherPropertiesPerSubstance;
  
                                                                                                      annotation (
      Icon(coordinateSystem(
            preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
          graphics={Text(
            extent={{-84,22},{92,64}},
            lineColor={0,0,255},
          textString="%name")}),
      Documentation(revisions="<html>
  <p>2009-2015 by Marek Matejak, Charles University, Prague, Czech Republic </p>
  </html>", info="<html>
  <h4>n = x &middot; n(solution) = &int; MolarFlow</h4>
  <p>where n is amount of the substance and x is mole fraction.</p>
  <p>The main class from &ldquo;Chemical&rdquo; package is called &quot;Substance&quot;. It has one chemical connector, where chemical potential and molar flow is presented. An amount of solute &quot;n&quot; is accumulated by molar flow inside an instance of this class. In the default setting the amount of solution &quot;n(solution)&quot; is set to 55.6 as amount of water in one liter, so in this setting the concentration of very diluted solution in pure water at &ldquo;mol/L&rdquo; has the same value as the amount of substance at &ldquo;mol&rdquo;. But in the advanced settings the default amount of solution can be changed by parameter or using solution port to connect with solution. The molar flow at the port can be also negative, which means that the solute leaves the Substance instance.&nbsp;</p>
  <p><br>The recalculation between mole fraction, molarity and molality can be written as follows:</p>
  <p>x = n/n(solution) = b * m(solvent)/n(solution) = c * V(solution)/n(solution)</p>
  <p>where m(solvent) is mass of solvent, V(solution) is volume of solution, b=n/m(solvent) is molality of the substance, c=n/V(solution) is molarity of the substance.</p>
  <p>If the amount of solution is selected to the number of total solution moles per one kilogram of solvent then the values of x will be the same as molality.</p>
  <p>If the amount of solution is selected to the number of total solution moles in one liter of solution then the values of x will be the same as molarity.</p>
  <p><br><br>Definition of electro-chemical potential:</p>
  <h4>u = u&deg; + R*T*ln(gamma*x) + z*F*v</h4>
  <h4>u&deg; = DfG = DfH - T * DfS</h4>
  <p>where</p>
  <p>x .. mole fraction of the substance in the solution</p>
  <p>T .. temperature in Kelvins</p>
  <p>v .. relative eletric potential of the solution</p>
  <p>z .. elementary charge of the substance (like -1 for electron, +2 for Ca^2+)</p>
  <p>R .. gas constant</p>
  <p>F .. Faraday constant</p>
  <p>gamma .. activity coefficient</p>
  <p>u&deg; .. chemical potential of pure substance</p>
  <p>DfG .. free Gibbs energy of formation of the substance</p>
  <p>DfH .. free enthalpy of formation of the substance</p>
  <p>DfS .. free entropy of formation of the substance </p>
  <p><br>Be carefull, DfS is not the same as absolute entropy of the substance S&deg; from III. thermodinamic law! It must be calculated from tabulated value of DfG(298.15 K) and DfH as DfS=(DfH - DfG)/298.15. </p>
  </html>"));
  end substanceb;

  model reaction "Chemical Reaction"
  
    parameter Integer nS=1 "Number of substrates types"
      annotation ( HideResult=true);
  
    parameter Modelica.SIunits.StoichiometricNumber s[nS]=ones(nS)
    "Stoichiometric reaction coefficient for substrates"
      annotation (HideResult=true);
  
    parameter Integer nP=1 "Number of products types"
      annotation ( HideResult=true);
  
    parameter Modelica.SIunits.StoichiometricNumber p[nP]=ones(nP)
    "Stoichiometric reaction coefficients for products"
      annotation (HideResult=true);
  
    Modelica.SIunits.MolarFlowRate rr(start=0) "Reaction molar flow rate";
    extends Chemical.Interfaces.ConditionalKinetics;
  
  Chemical.Interfaces.SubstancePorts_b substrates[nS] annotation (Placement(
        transformation(extent={{-110,-40},{-90,40}}), iconTransformation(
          extent={{-110,-40},{-90,40}})));
  Chemical.Interfaces.SubstancePorts_b products[nP] annotation (Placement(
        transformation(extent={{90,-40},{110,40}}), iconTransformation(extent=
           {{90,-40},{110,40}})));
  
    parameter Real kE(unit="mol/J")=0 "Kinetic turnover coefficient";
  
  protected
    Modelica.SIunits.ChemicalPotential du;
  equation
    //the main equation
    du = ((p * products.u) - (s * substrates.u));
    rr = kC * du * exp(-kE*abs(du));
  
    //reaction molar rates
    rr*s = -substrates.q;
    rr*p = products.q;
  
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
            100,100}}),   graphics={
          Rectangle(
            extent={{-100,-30},{100,30}},
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-100,-72},{100,-40}},
            lineColor={0,0,255},
          textString="%name"),
          Polygon(
            points={{-60,6},{-60,4},{54,4},{54,4},{18,14},{18,6},{-60,6}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{54,-8},{54,-6},{-60,-6},{-60,-6},{-24,-16},{-24,-8},{54,-8}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}),
      Documentation(revisions="<html>
  <p><i>2013-2015 by </i>Marek Matejak, Charles University, Prague, Czech Republic </p>
  </html>",     info="<html>
  <p><b>s<sub>1</sub>&middot;S<sub>1</sub> + .. + s<sub>nS</sub>&middot;S<sub>nS</sub> &lt;-&gt; p<sub>1</sub>&middot;P<sub>1</sub> + .. + p<sub>nP</sub>&middot;P<sub>nP</sub></b> </p>
  <p>By redefinition of stoichometry as v<sub>i</sub> = -s<sub>i</sub>, A<sub>i</sub> = S<sub>i</sub> for i=1..nS v<sub>i</sub> = p<sub>i-nS</sub>, A<sub>i</sub> = P<sub>i-nS</sub> for i=nS+1..nS+nP </p>
  <p>So the reaction can be written also as 0 = &sum; (v<sub>i</sub> &middot; A<sub>i</sub>) </p>
  <h4><span style=\"color:#008000\">Equilibrium equation</span></h4>
  <table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
  <td><p>K = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(S)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s) / <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>( a(P)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>s ) = <a href=\"modelica://ModelicaReference.Operators.'product()'\">product</a>(a(A)<a href=\"ModelicaReference.Operators.ElementaryOperators\">.^</a>v)&nbsp;</p></td>
  <td><p>dissociation constant</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>r</sub>G = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>G<sub>i</sub>) = &Delta;<sub>r</sub>H - T&middot;&Delta;<sub>r</sub>S = -R&middot;T&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(K) </p></td>
  <td><p>molar Gibb&apos;s energy of the reaction</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>r</sub>H = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>H<sub>i</sub>) </p></td>
  <td><p>molar enthalpy of the reaction</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>r</sub>S = &sum; (v<sub>i</sub> &middot; &Delta;<sub>f</sub>S<sub>i</sub>) = <a href=\"modelica://Modelica.Constants\">k</a>&middot;<a href=\"modelica://ModelicaReference.Operators.'log()'\">log</a>(&Delta;<sub>r</sub>&omega;) </p></td>
  <td><p>molar entropy of the reaction</p></td>
  </tr>
  </table>
  <h4><span style=\"color:#008000\">Notations</span></h4>
  <table cellspacing=\"2\" cellpadding=\"0\" border=\"0\"><tr>
  <td><p>A<sub>i</sub></p></td>
  <td><p>i-th substance</p></td>
  </tr>
  <tr>
  <td><p>v<sub>i</sub></p></td>
  <td><p>stochiometric coefficients of i-th substance</p></td>
  </tr>
  <tr>
  <td><p>K</p></td>
  <td><p>dissociation constant (activity based)</p></td>
  </tr>
  <tr>
  <td><p>a(A<sub>i</sub>)=f<sub>i</sub>*x<sub>i</sub></p></td>
  <td><p>activity of the substance A</p></td>
  </tr>
  <tr>
  <td><p>f<sub>i</sub></p></td>
  <td><p>activity coefficient of the substance A</p></td>
  </tr>
  <tr>
  <td><p>x<sub>i</sub></p></td>
  <td><p>mole fraction of the substance A</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>f</sub>H<sub>i</sub></p></td>
  <td><p>molar enthalpy of formation of i-th substance</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>f</sub>G<sub>i</sub></p></td>
  <td><p>molar Gibbs energy of formation of i-th substance</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>f</sub>S<sub>i</sub></p></td>
  <td><p>molar entropy of formation of i-th substance</p></td>
  </tr>
  <tr>
  <td><p>&Delta;<sub>r</sub>&omega;</p></td>
  <td><p>change of number of microstates of particles by reaction</p></td>
  </tr>
  <tr>
  <td></td>
  <td></td>
  </tr>
  </table>
  </html>"));
  end reaction;

  model esteva
    proyecto.solution solution1(ConstantTemperature = false)  annotation(
      Placement(visible = true, transformation(origin = {2, 0}, extent = {{-100, -100}, {100, 100}}, rotation = 0)));
    reaction reaction1 annotation(
      Placement(visible = true, transformation(origin = {-2, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  substance substance1 annotation(
      Placement(visible = true, transformation(origin = {-48, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  proyecto.substanceb substanceb1 annotation(
      Placement(visible = true, transformation(origin = {54, -2}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  equation
    connect(substanceb1.solution, solution1.solution) annotation(
      Line(points = {{60, -12}, {62, -12}, {62, -98}, {62, -98}}, color = {127, 127, 0}));
    connect(substance1.solution, solution1.solution) annotation(
      Line(points = {{-54, -8}, {-54, -8}, {-54, -60}, {62, -60}, {62, -98}, {62, -98}}, color = {127, 127, 0}));
    connect(substanceb1.port_a, reaction1.products[1]) annotation(
      Line(points = {{44, -2}, {8, -2}, {8, 0}, {8, 0}}, color = {158, 66, 200}));
    connect(substance1.port_a, reaction1.substrates[1]) annotation(
      Line(points = {{-38, 2}, {-12, 2}, {-12, 0}, {-12, 0}}, color = {158, 66, 200}));
  
  end esteva;











end proyecto;
