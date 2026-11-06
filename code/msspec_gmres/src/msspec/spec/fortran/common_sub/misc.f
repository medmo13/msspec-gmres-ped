      SUBROUTINE CLOSE_ALL_FILES()
          IMPLICIT NONE
          LOGICAL Ok
          INTEGER U

          DO U=7,200
              INQUIRE(UNIT=U, OPENED=Ok)
              IF (Ok) THEN
                  CLOSE(U)
              ENDIF
          ENDDO
      END SUBROUTINE CLOSE_ALL_FILES
