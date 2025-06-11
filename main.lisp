(defconstant +foot-to-metre+ 0.3048)

;; gas = (nitrogen helium oxygen)
(defparameter *nitrogen-index* 0)
(defparameter *oxygen-index* 1)
(defparameter *helium-index* 2)

;; common gas mixes
(defparameter *air*    '(.79 .00 .21))
(defparameter *eanx32* '(.68 .00 .32))
(defparameter *eanx36* '(.64 .00 .36))

(defparameter *model* 'zh-l16a)


(defun switch-model (model)
  ;; reinitialise everything -- can't be done mid-dive since number of compartments may differ
  (setf *model* model))


(defparameter *halftimes-matrix*
  '((zh-l16a . ((4.0 1.5)
                (8.0 3.0)
                (12.5 4.7)
                (18.5 7.0)
                (27.0 10.2)
                (38.3 14.5)
                (54.3 20.5)
                (77.0 29.1)
                (109.0 41.1)
                (146.0 55.1)
                (187.0 70.6)
                (239.0 90.2)
                (305.0 115.1)
                (390.0 147.2)
                (498.0 187.9)
                (635.0 239.6)))

    (zh-l16c . ((5.0    1.88)
                (8.0    3.02)
                (12.5   4.72)
                (18.5   6.99)
                (27.0   10.21)
                (38.3   14.48)
                (54.3   20.53)
                (77.0   29.11)
                (109.0  41.20)
                (146.0  55.19)
                (187.0  70.69)
                (239.0  90.34)
                (305.0  115.29)
                (390.0  147.42)
                (498.0  188.24)
                (635.0  240.03)))))


(defun get-model-halftimes (model)
  (rest (assoc model *halftimes-matrix*)))


(defun get-gas-halftimes (model gas-idx)
  (mapcar (lambda (compartment) (nth gas-idx compartment)) (get-model-halftimes model)))


(defun gas-mix-idx (gas-mix idx) (nth idx gas-mix))

(defun nitrogen-component (gas-mix) (gas-mix-idx gas-mix *nitrogen-index*))
(defun helium-component (gas-mix) (gas-mix-idx gas-mix *helium-index*))
(defun oxygen-component (gas-mix) (gas-mix-idx gas-mix *oxygen-index*))


(defun feet->metres (feet) (* feet +foot-to-metre+))
(defun metres->feet (metres) (/ metres +foot-to-metre+))


(defun depth->pressure (depth)
  "bar at metres"
  (1+ (/ depth 10)))


(defun pressure->depth (pressure)
  (* 10 (1- pressure)))


(defun gas-partial-pressure-at-pressure (gas-mix index pressure)
  (* pressure (nth index gas-mix)))


(defun gas-partial-pressure-at-depth (gas-mix index depth)
  (gas-partial-pressure-at-pressure
   gas-mix index (depth->pressure depth)))


(defun gas-compartment-halftime (inert-gas-index compartment model)
  (nth inert-gas-index (nth compartment (get-model-halftimes model))))


;; Pcomp = Pbegin + [ Pgas - Pbegin ] x [ 1 - 2 ^ ( - te / tht ) ]
;; where:
;; Pbegin = Inert gas pressure in the compartment before the exposure time ( ATM )
;; Pcomp = Inert gas pressure in the compartment after the exposure time ( ATM )
;; Pgas = Inert gas pressure in the mixture being breathed ( ATM )
;; te = Length of the exposure time ( minutes )
;; tht = Half time of the compartment ( minutes )

(defun tissue-pressure (pbegin pgas te tht)
  (+ pbegin
     (* (- pgas pbegin)
        (- 1 (expt 2 (/ (- te) tht))))))

(defun tissue-pressure-at-depth
    (compartment gas-mix gas-index start-pressure depth exposure-time model)
  "Partial pressure of gas in tissue compartment when at depth for exposure-time"
    (tissue-pressure
     start-pressure
     (gas-partial-pressure-at-depth gas-mix gas-index depth)
     exposure-time
     (gas-compartment-halftime gas-index compartment model)))


(defun tissue-pressures-at-depth
  (gas-mix gas-index initial-pressures depth exposure-time model)
  "Partial pressures of gas in all tissue compartments"

  (let ((compartments (loop for i from 0 repeat (length initial-pressures) collect i)))
    (mapcar (lambda (compartment pressure)
              (tissue-pressure-at-depth
               compartment
               gas-mix
               gas-index
               pressure
               depth
               exposure-time
               model))
            compartments
            initial-pressures)))


#| EXAMPLE:
    A diver descends from the surface to 100 feet on air and waits there
    ten minutes. The partial pressure of nitrogen in the breathing gas Pgas
    is 4 x 0.79 = 3.16 ATM. Let's pick a compartment, say number five.
    The nitrogen half-time for compartment five tht is 27 minutes. The
    nitrogen partial pressure in compartment five on the surface Pbegin is
    0.79 ATM, assuming the diver hasn't already been diving or subject to
    any altitude changes. The length of the exposure te is ten minutes.

    ANSWER: Pcomp = 1.33 bar
|#

(let* ((depth 30)
       (nitrogen-compartment-halftime 27)
       (exposure-time 10)

       (pbegin (nitrogen-component *air*))
       (pgas (gas-partial-pressure-at-depth *AIR* 0 depth)))

  (format t "~%Pcomp = ~A + (~A - ~A) x (1 - 2^(~A/~A)) = ~A~%"
          pbegin
          pgas pbegin
          (- exposure-time)
          nitrogen-compartment-halftime

          (tissue-pressure pbegin pgas exposure-time nitrogen-compartment-halftime)))


#| Simpler |#
(let ((pbegin (nitrogen-component *AIR*)))
  (print (tissue-pressure-at-depth 4 *AIR* *NITROGEN-INDEX* pbegin 30 10 'zh-l16a)))


#| Let's do it for each of the 16 tissue compartments, now under ZH-L16C |#
(let ((surface-load (loop repeat 16 collect (nitrogen-component *AIR*))))
  (print (tissue-pressures-at-depth *AIR* *NITROGEN-INDEX* surface-load 30 10 'zh-l16c)))



(defun halftime-a (tht)
  (* 2 (expt tht (- (/ 3)))))


(defun halftime-b (tht)
  (- 1.005 (expt tht (- (/ 2)))))


(defun a-b-parameters (tht)
  (values (halftime-a tht) (halftime-b tht)))


#| EXAMPLE:
For example, the half-time for compartment 5 is 27 minutes, so
a = 2 x ( 27 ^ -1/3 ) = 0.6667
b = 1.005 - ( 27^ -1/2 ) = 0.8125
|#

(let* ((compartment 4)
       (tht (gas-compartment-halftime *nitrogen-index* compartment 'zh-l16a)))

  (print (halftime-a tht))
  (print (halftime-b tht)))



#| Pambtol = ( Pcomp - a ) x b
where:
Pcomp = the inert gas pressure in the compartment ( ATM )
Pambtol = is the pressure you could drop to ( ATM )

EXAMPLE: Continuing the example above, we found that a exposure for ten
minutes to 4 ATM pressure ( 100 feet depth ), led to a nitrogen
pressure of 1.33 ATM in compartment 5. The a and b values for
compartment 5 were 0.6667 and 0.8125 respectively. Plugging these
into the above gives:

Pambtol = ( 1.33 - 0.6667 ) x 0.8125 = 0.54 ATM
|#

(defun compartment-decompression-pressure-ceiling (compartment pcomp model)
  "Pressure to which we can ascend given current compartment's inert gas pressure"
  (multiple-value-bind (a b)
      (a-b-parameters (gas-compartment-halftime *nitrogen-index* compartment model))

    (* b (- pcomp a))))


(defun decompression-ceiling (compartment-decompression-ceilings)
  (reduce #'max compartment-decompression-ceilings))


(defun no-stop-dive-p (compartment-decompression-ceilings)
  (< (decompression-ceiling compartment-decompression-ceilings) 1))


#| Compute the decompression ceiling in metres for compartment 5 when at 1.33bar |#
(let ((ceiling (pressure->depth (compartment-decompression-pressure-ceiling 4 1.33 'zh-l16a))))
  (print ceiling))


#| Compute the overall decompression ceiling for a given tissue load |#
(loop with exposure-time = 50
      with depth = 30
      with model = 'zh-l16a
      with start-pressure = (nitrogen-component *AIR*)

      initially (format t "~%Diving to 30m for 50 minutes~%")
      initially (format t "~%Tissue ~T Pressure (bar) ~T Ceiling (m)~%")
      for compartment from 0 upto 15

      for compartment-pressure = (tissue-pressure-at-depth
                                  compartment
                                  *AIR*
                                  *NITROGEN-INDEX*
                                  start-pressure
                                  depth
                                  exposure-time
                                  model)

      for compartment-ceiling = (compartment-decompression-pressure-ceiling
                                 compartment
                                 compartment-pressure
                                 model)

      do (format t "~3D ~10T ~4F ~30T ~4F~%"
                 compartment
                 compartment-pressure
                 (pressure->depth compartment-ceiling))

      collect compartment-ceiling into ceilings

      finally (format t "~%Decompression ceiling: ~4Fm (no stop dive? ~A)~%"
                      (pressure->depth (decompression-ceiling ceilings))
                      (if (no-stop-dive-p ceilings) "Y" "N")))


#| Multi-level dive
    Starting at the surface, we then descend to 30m for 10 minutes, followed by 40m for 5 minutes
|#

(let* ((compartments (loop for i from 0 upto 15 collect i))
       (tissue-load (mapcar (lambda (compartment)
                              (declare (ignore compartment))
                              (nitrogen-component *AIR*))
                            compartments)))

  (format t "~%Initial load: ~A~%" tissue-load)

  (setf tissue-load (tissue-pressures-at-depth
                     *AIR* *NITROGEN-INDEX* tissue-load 30 10 'zh-l16a))
  (format t "~%Descend to 10m for 5 minutes: ~A~%" tissue-load)

  (setf tissue-load (tissue-pressures-at-depth
                     *AIR* *NITROGEN-INDEX* tissue-load 40 5 'zh-l16a))
  (format t "~%Descend to 10m for 5 minutes: ~A~%" tissue-load))


#| Observe that tissue loading is less if we just descend straight to 40m for 5 minutes |#
(let* ((compartments (loop for i from 0 upto 15 collect i))
       (tissue-load (mapcar (lambda (compartment)
                              (declare (ignore compartment))
                              (nitrogen-component *AIR*))
                            compartments)))

  (format t "~%Initial load: ~A~%" tissue-load)

  (setf tissue-load (tissue-pressures-at-depth
                     *AIR* *NITROGEN-INDEX* tissue-load 40 5 'zh-l16a))
  (format t "~%Descend to 40m for 5 minutes: ~A~%" tissue-load))
