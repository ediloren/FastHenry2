' FastHenry2 Automation example
' Enrico Di Lorenzo, 2004/05/07
'
' Launch from a DOS shell, typing 'cscript <filename>',
' where <filename> is the name of this file
' (can also double click on the file, but you only get
' message boxes in this case)

Dim FastHenry2
' Create FastHenry2 object
Set FastHenry2 = CreateObject("FastHenry2.Document")

' Extract script path from ScriptFullName Property
pathPos = InstrRev(Wscript.ScriptFullName, Wscript.ScriptName)
path = left(Wscript.ScriptFullName, pathPos-1)

For  i= 1 to 3
  ' Input files are already existing in this simple example;
  ' could also generate them on the fly and stop loop when
  ' desired refinement accuracy is reached (e.g. <1% difference
  ' between FastHenry2 results) or optimization result is achieved
  
  ' Try to run FastHenry2
  ' Remark: the run path must be surrounded by quotas '"' to support
  ' also paths containing spaces (quotas in VisualBasic are escaped by
  ' doubling the symbol, i.e., "" )
  couldRun = FastHenry2.Run("""" + path + "coils"+CStr(i)+".inp""")
  ' Wait for end of operation, using polling; could also use callback function
  ' see Windows Script Components help
  Do While FastHenry2.IsRunning = True
    Wscript.Sleep 500
  Loop
  ' retrieve inductance matrix
  inductance = FastHenry2.GetInductance()
  WScript.Echo "Coils"+CStr(i)+" mutual inductance is " + CStr(inductance(0, 0, 1))
Next

' Quit FastHenry2
FastHenry2.Quit
' Destroy FastHenry2 object
Set FastHenry2 = Nothing
