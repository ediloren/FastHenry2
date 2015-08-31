VERSION 1.0 CLASS
BEGIN
  MultiUse = 0   'False
  Persistable = 0  'NotPersistable
  DataBindingBehavior = 0  'vbNone
  DataSourceBehavior  = 0  'vbNone
  MTSTransactionMode  = 0  'NotAnMTSObject
END
Attribute VB_Name = "Callback"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = True
Attribute VB_PredeclaredId = False
Attribute VB_Exposed = True
Public Sub FastHenryEndCallback()
    ' signal end of FastHenry2 execution
    mainForm.TextStatusWindow.Text = "FastHenry2 STOPPED"
End Sub

Public Sub FastHenryLogCallback(logString As String, color As Long)
    
    ' log messages from FastHenry2
    
    ' format carriage returns & line feeds
    logString = Replace(logString, Chr$(10), Chr$(13) + Chr$(10))
    
    ' trim output text if too long
    If Len(mainForm.TextLogWindow.Text) + Len(logString) > 15000 Then
        mainForm.TextLogWindow.Text = Right$(mainForm.TextLogWindow.Text, Len(mainForm.TextLogWindow.Text) - Len(logString) - 5000)
    End If
       
    ' output text
    mainForm.TextLogWindow.SelLength = 0
    mainForm.TextLogWindow.SelStart = Len(mainForm.TextLogWindow.Text)
    mainForm.TextLogWindow.SelColor = color
    mainForm.TextLogWindow.SelText = logString
    
    'mainForm.TextLogWindow.Text = mainForm.TextLogWindow.Text + Replace(logString, Chr$(10), Chr$(13) + Chr$(10))
End Sub
