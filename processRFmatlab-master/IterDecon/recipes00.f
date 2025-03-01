      SUBROUTINE CORREL(DATA1,DATA2,N,ANS)
      PARAMETER(NMAX=8192)
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX FFT(NMAX),ANS(N)
      CALL TWOFFT(DATA1,DATA2,FFT,ANS,N)
      NO2=FLOAT(N)/2.0
      DO 11 I=1,N/2+1
        ANS(I)=FFT(I)*CONJG(ANS(I))/NO2
11    CONTINUE
      ANS(1)=CMPLX(REAL(ANS(1)),REAL(ANS(N/2+1)))
      CALL REALFT(ANS,N/2,-1)
      RETURN
      END
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
      ELSE
        C2=0.5
        THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 11 I=2,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
      SUBROUTINE TWOFFT(DATA1,DATA2,FFT1,FFT2,N)
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX FFT1(N),FFT2(N),H1,H2,C1,C2
      C1=CMPLX(0.5,0.0)
      C2=CMPLX(0.0,-0.5)
      DO 11 J=1,N
        FFT1(J)=CMPLX(DATA1(J),DATA2(J))
11    CONTINUE
      CALL FOUR1(FFT1,N,1)
      FFT2(1)=CMPLX(AIMAG(FFT1(1)),0.0)
      FFT1(1)=CMPLX(REAL(FFT1(1)),0.0)
      N2=N+2
      DO 12 J=2,N/2+1
        H1=C1*(FFT1(J)+CONJG(FFT1(N2-J)))
        H2=C2*(FFT1(J)-CONJG(FFT1(N2-J)))
        FFT1(J)=H1
        FFT1(N2-J)=CONJG(H1)
        FFT2(J)=H2
        FFT2(N2-J)=CONJG(H2)
12    CONTINUE
      RETURN
      END
