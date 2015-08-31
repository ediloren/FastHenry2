VERSION 5.00
Object = "{3B7C8863-D78F-101B-B9B5-04021C009402}#1.2#0"; "RICHTX32.OCX"
Begin VB.Form FastHenryClient 
   Caption         =   "FastHenry2 Test Driver"
   ClientHeight    =   6345
   ClientLeft      =   1635
   ClientTop       =   2190
   ClientWidth     =   6735
   BeginProperty Font 
      Name            =   "MS Sans Serif"
      Size            =   8.25
      Charset         =   0
      Weight          =   700
      Underline       =   0   'False
      Italic          =   0   'False
      Strikethrough   =   0   'False
   EndProperty
   ForeColor       =   &H80000008&
   LinkTopic       =   "Form1"
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   6345
   ScaleWidth      =   6735
   Begin VB.ComboBox ComboSelectMatrix 
      Height          =   315
      ItemData        =   "fhdriv.frx":0000
      Left            =   5160
      List            =   "fhdriv.frx":000A
      Style           =   2  'Dropdown List
      TabIndex        =   17
      Top             =   2880
      Width           =   1335
   End
   Begin RichTextLib.RichTextBox TextLogWindow 
      Height          =   1335
      Left            =   240
      TabIndex        =   16
      Top             =   4320
      Width           =   6255
      _ExtentX        =   11033
      _ExtentY        =   2355
      _Version        =   393217
      Enabled         =   -1  'True
      ScrollBars      =   3
      TextRTF         =   $"fhdriv.frx":0026
   End
   Begin VB.CheckBox CheckLogCallback 
      Caption         =   "Log Callback"
      Height          =   375
      Left            =   4560
      TabIndex        =   14
      Top             =   3360
      Width           =   1695
   End
   Begin VB.CheckBox CheckEndCallback 
      Caption         =   "End Callback"
      Height          =   375
      Left            =   4560
      TabIndex        =   13
      Top             =   3720
      Width           =   1695
   End
   Begin VB.CommandButton CommandStop 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "&Stop FH"
      Height          =   375
      Left            =   4560
      TabIndex        =   11
      Top             =   1920
      Width           =   1935
   End
   Begin VB.CommandButton CommandShow 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "Show Window"
      Height          =   375
      Left            =   4560
      TabIndex        =   9
      Top             =   480
      Width           =   1935
   End
   Begin VB.TextBox TextCmdLine 
      Appearance      =   0  'Flat
      Height          =   375
      Left            =   240
      ScrollBars      =   1  'Horizontal
      TabIndex        =   8
      Text            =   "C:\temp\fastmodel\fasthenry-3.0\examples\broken.inp"
      Top             =   3360
      Width           =   4095
   End
   Begin VB.CommandButton CommandRun 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "&Run FH"
      Height          =   375
      Left            =   4560
      TabIndex        =   7
      Top             =   960
      Width           =   1935
   End
   Begin VB.TextBox TextStatusWindow 
      Appearance      =   0  'Flat
      Height          =   375
      Left            =   240
      TabIndex        =   6
      Top             =   2640
      Width           =   4095
   End
   Begin VB.CommandButton CloseFH 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "Quit FH"
      Height          =   375
      Left            =   4560
      TabIndex        =   3
      Top             =   2400
      Width           =   1935
   End
   Begin VB.CommandButton GetMatrix 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "Get"
      Height          =   375
      Left            =   4560
      TabIndex        =   4
      Top             =   2880
      Width           =   495
   End
   Begin VB.CommandButton IsRunning 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "&IsRunning"
      Height          =   375
      Left            =   4560
      TabIndex        =   2
      Top             =   1440
      Width           =   1935
   End
   Begin VB.TextBox TextResult 
      Appearance      =   0  'Flat
      BeginProperty Font 
         Name            =   "Courier New"
         Size            =   8.25
         Charset         =   0
         Weight          =   400
         Underline       =   0   'False
         Italic          =   0   'False
         Strikethrough   =   0   'False
      EndProperty
      Height          =   1815
      Left            =   240
      MultiLine       =   -1  'True
      ScrollBars      =   3  'Both
      TabIndex        =   1
      Top             =   480
      Width           =   4095
   End
   Begin VB.CommandButton Close 
      Appearance      =   0  'Flat
      BackColor       =   &H80000005&
      Caption         =   "&Close"
      Height          =   375
      Left            =   4800
      TabIndex        =   5
      Top             =   5880
      Width           =   1695
   End
   Begin VB.Label Label4 
      Caption         =   "Log Window"
      Height          =   255
      Left            =   240
      TabIndex        =   15
      Top             =   3960
      Width           =   1335
   End
   Begin VB.Label Label2 
      Caption         =   "Command Line"
      Height          =   255
      Left            =   240
      TabIndex        =   12
      Top             =   3120
      Width           =   1335
   End
   Begin VB.Label Label1 
      Caption         =   "Status Window"
      Height          =   255
      Left            =   240
      TabIndex        =   10
      Top             =   2400
      Width           =   1335
   End
   Begin VB.Label Label3 
      Caption         =   "Resistance / Inductance"
      Height          =   255
      Left            =   240
      TabIndex        =   0
      Top             =   240
      Width           =   2535
   End
End
Attribute VB_Name = "FastHenryClient"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Dim FastHenry2 As Object
Dim callbackObj As New Callback
Dim NL As String

Private Sub CheckEndCallback_Click()
    ' If checked
    If CheckEndCallback.Value = 1 Then
        ' Set a callback function to be called upon FastHenry2 end of simulation
        vResult = FastHenry2.SetEndCallback(callbackObj, "FastHenryEndCallback")
    Else
        ' Set no callback
        vResult = FastHenry2.SetEndCallback(Nothing, "FastHenryEndCallback")
    End If
End Sub

Private Sub CheckLogCallback_Click()
    ' If checked
    If CheckLogCallback.Value = 1 Then
        ' Set a callback function to be called upon FastHenry2 end of simulation
        vResult = FastHenry2.SetLogCallback(callbackObj, "FastHenryLogCallback")
    Else
        ' Set no callback
        vResult = FastHenry2.SetLogCallback(Nothing, "FastHenryLogCallback")
    End If
End Sub

Private Sub Close_Click()
    FastHenry2.Quit
    End
End Sub

Private Sub CloseFH_Click()
    ' Try to quit FastHenry2
    couldQuit = FastHenry2.Quit
    If couldQuit Then
        TextStatusWindow.Text = "FastHenry2 QUITTED"
    Else
        TextStatusWindow.Text = "FastHenry2 NOT quitted"
    End If
End Sub

Private Sub CommandRun_Click()
    ' Try to run FastHenry2
    couldRun = FastHenry2.Run(TextCmdLine.Text)
    If couldRun Then
        TextStatusWindow.Text = "FastHenry2 STARTED"
    Else
        TextStatusWindow.Text = "FastHenry2 NOT started"
    End If

End Sub

Private Sub CommandShow_Click()
    ' Show FastHenry2 window
    FastHenry2.ShowWindow
End Sub

Private Sub CommandStop_Click()
    ' Force stop of simulation
    FastHenry2.Stop
End Sub

Private Sub Form_Load()
    ' Create FastHenry2 object
    Set FastHenry2 = CreateObject("FastHenry2.Document")
    ComboSelectMatrix.ListIndex = 0
    NL = Chr$(13) + Chr$(10)
End Sub

Private Sub Form_Unload(Cancel As Integer)
    ' Destroy FastHenry2 object
    Set FastHenry2 = Nothing
End Sub

Private Sub GetMatrix_Click()
    TextResult.Text = ""
    ' retrieve port names
    rowPortNames = FastHenry2.GetRowPortNames()
    colPortNames = FastHenry2.GetColPortNames()
    ' get number or rows, cols
    rowDim = UBound(rowPortNames, 1)
    colDim = UBound(colPortNames, 1)
    ' retrieve frequencies array
    frequencies = FastHenry2.GetFrequencies()
    freqDim = UBound(frequencies, 1)
   
    ' if there are as many rows as columns, all ok
    If rowDim = colDim Then
        ' retrieve matrix
        If ComboSelectMatrix.ListIndex = 0 Then
            TextResult.Text = "Resistance matrix (Ohm)" + NL + NL
            matrix = FastHenry2.GetResistance()
        Else
            TextResult.Text = "Inductance matrix (Henry)" + NL + NL
            matrix = FastHenry2.GetInductance()
        End If
        ' output matrix values
        For k = 0 To freqDim
            TextResult.Text = TextResult.Text + "Frequency: " + Format(frequencies(k), "#.00E+00") + NL + NL
            For i = 0 To rowDim
                TextResult.Text = TextResult.Text + Format(Left$(rowPortNames(i), 12), "!@@@@@@@@@@@@") + " "
                For j = 0 To rowDim
                    TextResult.Text = TextResult.Text + Format(matrix(k, i, j), "\ #.000E+00;-#.000E+00") + " "
                Next j
                TextResult.Text = TextResult.Text + NL
            Next i
            TextResult.Text = TextResult.Text + NL
        Next k
    ' if number of columns <> number of rows, option '-x' has been used
    ' in FastHenry2. In this case, FastHenry2 computes an ADMITTANCE
    ' rectangular matrix, see FastHenry2 user's manual.
    Else
        ' retrieve matrix
        If ComboSelectMatrix.ListIndex = 0 Then
            TextResult.Text = "Admittance matrix, real part" + NL + NL
            matrix = FastHenry2.GetResistance()
        Else
            TextResult.Text = "Admittance matrix, imaginary part" + NL + NL
            matrix = FastHenry2.GetInductance()
        End If
        ' output matrix values
        For k = 0 To freqDim
            TextResult.Text = TextResult.Text + "Frequency: " + Format(frequencies(k), "#.00E+00") + NL + NL
            ' column names
            TextResult.Text = TextResult.Text + "             "
            For j = 0 To colDim
                TextResult.Text = TextResult.Text + Format(Left$(colPortNames(j), 9), "!@@@@@@@@@") + " "
            Next j
            TextResult.Text = TextResult.Text + NL
            For i = 0 To rowDim
                TextResult.Text = TextResult.Text + Format(Left$(rowPortNames(i), 12), "!@@@@@@@@@@@@") + " "
                For j = 0 To colDim
                    TextResult.Text = TextResult.Text + Format(matrix(k, i, j), "\ #.000E+00;-#.000E+00") + " "
                Next j
                TextResult.Text = TextResult.Text + NL
            Next i
            TextResult.Text = TextResult.Text + NL
        Next k
    End If
End Sub

Private Sub IsRunning_Click()
    ' Test if FastHenry2 is presently running a simulation
    If FastHenry2.IsRunning Then
        TextStatusWindow.Text = "FastHenry2 RUNNING"
    Else
        TextStatusWindow.Text = "FastHenry2 NOT running"
    End If
End Sub

