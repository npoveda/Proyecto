package proyecto
  model Solution "Chemical solution as homogenous mixture of the substances"
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
  end Solution;

  model Solucion "Chemical solution as homogenous mixture of the substances"
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
  end Solucion;





end proyecto;
