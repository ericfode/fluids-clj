(ns fluids.core
  (:require
   [fluids.gen :as g]

   [clojure.spec.alpha :as s]
   [clojure.spec.test.alpha :as st]
   [clojure.spec.gen.alpha :as gen]
   [uncomplicate.neanderthal.core :as nc]
   [uncomplicate.fluokitten.core :as fc]
   [uncomplicate.neanderthal.real :as nr]
   [uncomplicate.neanderthal.native :as nn]
   [uncomplicate.neanderthal.math :as nm]

   ))



(comment
"""
!*******************************************************************************
!*
!* Educationally-Designed Unstructured 1D (EDU1D) Code
!*
!*  --- EDU1D oned_first_order_diffusion
!*
!*
!* One-dimensional first-order Diffusion Scheme.
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 1 (05-26-2018).
!*
!* -05-26-18: A typo corrected for Lr: two**pi -> two*pi
!*  thanks to Sinath at University of Tokyo.
!*
!* This F90 program was written and made available for download
!* for an educational purpose. Comments are welcome.
!*
!* This file may be updated in future.
!*
!* Katate Masatsuka, September 2017. http://www.cfdbooks.com
!*
!*-----------------------------------------------------------------------------
!*
!*  Problem: u_tau = nu*u_{xx} + nu*pi^2*sin(pi*x), u(0)=u(1)=0,
!*
!*           which is solved in the first-order system form:
!*
!*           u_tau = nu*p_x + nu*pi^2*sin(pi*x)
!*           p_tau = (u_x - p)/Tr
!*
!*           The two formulations have the same steady solution.
!*
!* Note: The first-order system is hyperbolic in the pseudo time, tau.
!*       The eigenvalues (wave speeds) are + nu/Lr, and - nu/Lr.
!*
!* Note: The first-order system is equivalent to the original diffusion equation
!*       in the steady state. The idea is to integrate the first-order system in
!*       time, instead of the original diffusion equation. Thereby, the time step
!*       is O(h), not O(h^2), and accurate solution gradient can be computed
!*       simultaneously. See Nishikawa-JCP2014 for details.
!*
!* Note: Two schemes are implemented:
!*
!*       (1)First-order diffusion scheme based on the hyperbolic form.
!*          - Pseudo time step O(h)
!*          - Accuracy: error(u,p) = O(h)
!*          - Weak boundary condition.
!*          - This is Scheme-I in Ref.[1].
!*          - You can implement Scheme-II in Ref.[1], which is a better scheme.
!*
!*       (2)Conventional finite-difference scheme.
!*          - Pseudo time step O(h^2)
!*          - Accuracy: error(u  ) = O(h), error(p)=O(1).
!*          - Strong boundary condition.
!*
!* Note: Gradients computed from 1st-order solution are typically inconsistent
!*       on irregular grids. First-order scheme constructed by the hyperbolic
!*       method gives 1st-order accurate gradients. This is one of the advantages
!*       of the hyperbolic method, which leads to a usable 1st-order scheme
!*       for the Navier-Stokes equations. See discussions for the first-order
!*       accurate Navier-Stokes scheme in AIAA Paper-2091:
!*       http://hiroakinishikawa.com/My_papers/nishikawa_aiaa2014-2091.pdf
!*
!*
!* [1] H. Nishikawa, First-, second-, and third-order finite-volume schemes
!*     for diffusion, Journal of Computational Physics, 256, pp. 791-805, 2014
!*     https://doi.org/10.1016/j.jcp.2013.09.024
!* http://hiroakinishikawa.com/My_papers/nishikawa_jcp2014v256pp791-805_preprint.pdf
!*
!*******************************************************************************
 program oned_first_order_diffusion_scheme
""")

(comment """
 implicit none

  integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: zero=0.0_dp, one=1.0_dp, two=2.0_dp, half=0.5_dp
 real(dp), parameter :: pi=3.141592653589793238_dp

    real(dp), dimension(:), allocatable ::      x          ! Nodal coordinates
    real(dp), dimension(:), allocatable ::    vol          ! Dual volume at nodes
    real(dp), dimension(:), allocatable ::      u,      p  ! Solution & Gradient
    real(dp), dimension(:), allocatable :: uexact, pexact  ! Exact solutions
    real(dp), dimension(:), allocatable ::  res_u,  res_p  ! Nodal Residuals

    real(dp) :: nu          ! Diffusion coefficient (constant)
    real(dp) :: h           ! Mesh spacing for an uniform mesh.
    real(dp) :: hmin, hmax  ! Min/Max mesh spacing for a perturbed irregular mesh.
    real(dp) :: heff        ! Effective (averaged) mesh spacing.
    real(dp) :: rn          ! Random number
    real(dp) :: Lr, Tr      ! Length scale and relaxation time.
    real(dp) :: dtau        ! Pseudo time step.
    real(dp) :: res_max     ! Maximum residual to check the convergence.
    integer  :: nnodes      ! Total number of nodes
    integer  :: scheme_type ! = 1 : First-order diffusion scheme.
                            ! = 2 : Conventional finite-difference scheme
    integer  :: k, j

    real(dp) ::   f(2)      ! Numerical flux
    real(dp) ::   uL, uR    ! Left and right solutions at interface.
    real(dp) ::   pL, pR    ! Left and right solution gradients at interface.
    integer  ::   iu, ip

  iu = 1
  ip = 2
""")

(s/def ::dv  nc/vctr?)

(s/def ::x ::dv)

(s/def ::vol ::dv)
(s/def ::u ::dv)
(s/def ::p ::dv)
(s/def ::uexact ::dv)
(s/def ::pexact ::dv)
(s/def ::res_u ::dv)
(s/def ::res_p ::dv)

(s/def ::nu double?)
(s/def ::h double?)
(s/def ::hmin double?)
(s/def ::hmax double?)
(s/def ::rn double?)
(s/def ::Lr double?)
(s/def ::Tr double?)
(s/def ::dtau double?)
(s/def ::res_max double?)
(s/def ::nnodes int?)
(s/def ::scheme #{::hyperbolic ::euler})

(s/def ::state (s/keys :req [::x ::vol ::u ::p ::uexact ::pexact ::res_u ::res_p ::nu ::hmin ::hmax ::Lr ::Tr ::dtau ::res_max ::nnodes ::scheme]))


(comment """
!--------------------------------------------------------------------------------
! Initialization (values have no meaning; they will be overwritten later.)

    Lr = one
    Tr = one
  dtau = one
!--------------------------------------------------------------------------------
! Diffusion coefficient

  nu = one

!--------------------------------------------------------------------------------
! Input: Select the scheme.

  do j = 1, 100

  write(*,*) " Type of scheme = ?"
  write(*,*) 
  write(*,*) "   1 -> First-order scheme solving the hyperbolic system"
  write(*,*) "   2 -> Conventional finite-difference scheme"
  read(*,*) scheme_type

  if (scheme_type /= 1 .and. scheme_type /= 2) then
   write(*,*) " Wrong number. Enter 1 or 2."
  else
   exit
  endif

  end do

!--------------------------------------------------------------------------------
! Input: Set the number of nodes.

  write(*,*) " The number of nodes = ?"
  read(*,*) nnodes

""")

(defn init-state [^double size ]
  (let [  Lr  1.0
          Tr  1.0
        dtau  1.0
          nu  1.0])


  size
  )



(comment """
!--------------------------------------------------------------------------------
! Allocate arrays.

   allocate(                     x(nnodes)) ! Nodal coordinates.
   allocate(                   vol(nnodes)) ! Dual volume
   allocate(     u(nnodes),      p(nnodes)) ! Solution & Gradient
   allocate(uexact(nnodes), pexact(nnodes)) ! Exact solutions
   allocate( res_u(nnodes),  res_p(nnodes)) ! Nodal Residuals

!--------------------------------------------------------------------------------
! Set up                           h
!          Grid ->  o-------o-----o-----o------o---o-----o---o
!                  j=1            j    j+1                 j=nnodes
""")
(s/fdef nudge
  :args (s/cat :x (s/with-gen
                    (s/and ::dv #(= (nc/dim %) 4))
                    #(g/ordered-dv-gen)))
  :ret (s/and ::dv #(= (nc/dim %) 4))
  :fn (fn [{:keys [args ret]}]
        (and
         ;;The neighborhood should still be sorted
         (< (nr/entry ret 0) (nr/entry ret 1))
         (< (nr/entry ret 1) (nr/entry ret 2)))))

(defn nudge [x]
  "given a neighborhood of points in a vector and a rand [1 2 3 r]
   rand is normalized to the size of the neighborhood
   nudges the vector [1 2+rand 3 r]

(nn/dge 3 10 (nn/dv (flatten (repeat 10 [1 2 3]))))
"
  (let [jm1  (nr/entry x 0)
        j    (nr/entry x 1)
        jp1  (nr/entry x 2)
        r    (nr/entry x 3)
        l    (* 0.5  (+ jp1 jm1))
        m    (* 0.25 (- r 0.5))
        r2    (* 0.5 m (- jp1 jm1))]
    (nr/entry! x 1 (+ l r2))))

(defn approx-uniform-step? [tolerance coll]
  (let [diffs
        (->> coll
             (partition 2 1)
             (map (partial apply -)))
        [max-diff min-diff] ((juxt max min) diffs)]
    (or (= max-diff min-diff)
        (< (- max-diff min-diff) tolerance))))


(s/fdef normalized-points
  :args (s/cat :size  (s/int-in 2 10000))
  :ret  (s/and (s/coll-of double?)
               #(= 0.0 (first %))
               #(= 1.0 (last %))
               #(approx-uniform-step? 0.0001 %))
  :fn   (fn [{:keys [args ret]}]
          (let [size (:size args)]
            (= (count ret) size))))

(defn normalized-points [size]
  (let [h (/ 1 (dec size))]
    (into []
          (fc/fmap
           (fn [x] (double (* h (- x 1))))
           (range 1 (inc size))))))

(defn points-gen  [&{:keys [min-size max-size]
                     :or {min-size 2 max-size 1000}}]
  (gen/fmap
   (fn [size]
     [(normalized-points size)
      (flatten [0.0  (repeatedly (- size 2) rand) 0.0])])
   (s/gen (s/int-in min-size max-size))))

(s/fdef offset-matrix
  :args (s/with-gen
          (s/cat :points (s/coll-of double?)
                 :rands (s/coll-of double?))
          #(points-gen))
  :ret  nc/matrix?
  :fn   (fn [{:keys [args ret]}]
          (let [points-size  (count (:points args))
                rands-size   (count (:rangs args))
                ret-size     (nc/dim ret)]
            (= points-size rands-size)
            (= (nc/row ret 1) (nn/dv (:points args)))
            (= (nc/row ret 3) (nn/dv (:rands args)))
            (= ret-size (* 4 points-size)))))

(defn offset-matrix [points rands]
  (let [size (count points)]
    (nn/dge 4 size (interleave (cons 0.0 (drop-last points))
                                points
                                (conj (rest points) 1.0)
                                rands))))

(defn irregular-grid [size]
  "Returns a normailized irregular grid"
  (let [points         (normalized-points size)
        rands          (flatten [0  (repeatedly (- size 2) rand) 0])
        matrix         (offset-matrix points rands)
        preturb-matrix (nc/submatrix matrix  0 1  4  (dec size))
        preturbed      (fc/fmap nudge (nc/cols preturb-matrix))]
    (nn/dv
     (cons 0.0
           (fc/fmap #(nr/entry % 1)
                    preturbed)))))


(irregular-grid 5)

(defn gradients [state]
"
  !Compute and store the exact solution and the gradient.

  do j = 1, nnodes
  uexact(j) = sin(pi*x(j))     ! Exact solution at xj
  pexact(j) = pi*cos(pi*x(j))  ! Exact gradient at xj
  u(j) = zero             ! Initial solution
  p(j) = zero             ! Initial gradient
  end do
 "
  (let [singles
        (map
         #(hash-map
           :uexact (nm/sin (* nm/pi %))
           :pexact (* nm/pi (nm/cos (* nm/pi %)))
           :v %
           :u  0
           :p  0)
         state)]

    (reduce (fn [old new]
              (-> old
                  (update :uexact conj (:uexact new))
                  (update :u conj (:u new))
                  (update :v conj (:v new))
                  (update :pexact conj (:pexact new))
                  (update :p conj (:p new))))
            {:uexact []
             :pexact []
             :v []
             :u []
             :p []}
            singles)))




(defn bands
  [state]
  [(nc/subvector state 1 (dec (nc/dim state)))
   (nc/subvector state 0 (dec (nc/dim state)))])


(defn dual-volumes [state]
  "
 !Compute and store the dual volume around each node.

   do j = 2, nnodes-1
    vol(j) = half*( x(j+1)-x(j-1) )
   end do
    vol(1     ) = half*( x(2)-x(1) )
    vol(nnodes) = half*( x(nnodes)-x(nnodes-1) ) "
  (let [vols
        (apply mapv
               (fn [xp1 xm1]
                 (* 1/2 (- xp1 xm1)))
               (bands state))]
    (conj vols (last vols))))


(def t (irregular-grid 10))


(bands (irregular-grid 10))

(dual-volumes (irregular-grid 10))

(defn mesh-spacing
  [volumes count]
  "
 !Compute the effective mesh spacing: h = L1(vol).

    heff = zero
   do j = 1, nnodes
    heff = heff + vol(j)
   end do
    heff = heff/real(nnodes,dp)

 !Compute the minimum and maximum mesh spacings.

    hmin = minval(vol)
    hmax = maxval(vol)
"
  (let [vs (nn/dv volumes)]
    {:heff
     (/
      (reduce
       (fn [heff vol]
         (+ heff vol)) 0 volumes)
      count)
     :hmin (nth volumes (nc/imin vs))
     :hmax (nth volumes (nc/imax vs))}))

(mesh-spacing (dual-volumes (irregular-grid 10))

              10)

(defn upwind-baby
  [x vols
   {:keys [hmin hmax heff] :as spacing}
   {:keys [u p uexact pexact] :as gradients}]


  (println "X")
  (clojure.pprint/pprint x)
  
  (println "Vols")
  (clojure.pprint/pprint vols)

  (println "u")
  (clojure.pprint/pprint u)

  (println "p")
  (clojure.pprint/pprint p)

  (println "grad")
  (clojure.pprint/pprint gradients)
  (let [iu 0
        ip 1
        nu 1.0
        u (nn/dv u)
        p (nn/dv p)
        Lr (/ 1.0 (* 2.0 nm/pi))
        Tr (/ (* Lr Lr ) nu)
        vols (nn/dv vols)
        nnodes (nc/dim x)
        dtau (* 0.99 (/ hmin (/ nu Lr)))
        res_u (nn/dv (take nnodes (repeat 0.0)))
        res_p (nn/dv (take nnodes (repeat 0.0)))
        f (nn/dv [0 0])]
    (loop [pesudotime (range 10000000)
           x x
           u u
           p p
           f f
           res_u res_u
           res_p res_p]
      ;; interior_face_loop : do j = 1, nnodes-1 ! j-th face between j and j+1.
      ;; ! Left and right states.
      ;; uL = u(j)
      ;; pL = p(j)

      ;; uR = u(j+1)
      ;; pR = p(j+1)

      ;; ! Compute the numerical flux (upwind flux)

      ;; f(iu) = half*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) )
      ;; f(ip) = half*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) )

      ;; ! Add to the left and subtract from the right for Res = f_{j+1/2}-f_{j=1/2}.

      ;; res_u(j  ) = res_u(j  ) + f(iu)
      ;; res_p(j  ) = res_p(j  ) + f(ip)

      ;; res_u(j+1) = res_u(j+1) - f(iu)
      ;; res_p(j+1) = res_p(j+1) - f(ip)

      ;; end do interior_face_loop

      (do
        (for [index (range nnodes)]
          (let [uL (nth u index)
                pL (nth p index)
                uR (nth u (inc index))
                pR (nth p (inc index))
                newiu (* 1/2 (+ (* nu pR) (* nu pL)
                                (/ nu (* Lr (- uR uL)))))
                newip (* 1/2 (+       uR         uL
                                      (/ nu (* Lr (- pR pL)))))

                new_res_u_l  (+ (nc/entry res_u index) newiu)
                new_res_p_l  (+ (nc/entry res_p index) newip)
                new_res_u_r  (- (nc/entry res_u (inc index)) newiu)
                new_res_p_r  (- (nc/entry res_u (inc index)) newip)
                ]
            (do
              (clojure.pprint/pprint u)
              (clojure.pprint/pprint p)
              (clojure.pprint/pprint f)
              (clojure.pprint/pprint res_u)
              (clojure.pprint/pprint res_p))
            (nc/entry! f iu newiu)
            (nc/entry! f ip newip)
            (nc/entry! res_u index new_res_u_l)
            (nc/entry! res_p index new_res_p_l)
            (nc/entry! res_u (inc index) new_res_u_r)
            (nc/entry! res_p (inc index) new_res_p_r)))

        ;; !------------------------------------------------------------
        ;; ! Fluxes through the domain boundaries to close the residuals.
        ;; ! Note: Weak boundary condition is applied.
        ;; !------------------------------------------------------------

        ;;   !-------------------------------------------
        ;;   ! Left boundary
        ;;   !-------------------------------------------

        ;;       j = 1

        ;;      uL = 0.0  !<- Boundary condition: u(0)=0.
        ;;      pL = p(j) !<- No boundary condition: copy from the right.

        ;;      uR = u(j)
        ;;      pR = p(j)

        ;;     !Compute the numerical flux (upwind flux)

        ;;       f(iu) = half*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) )
        ;;       f(ip) = half*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) )

        ;;     !Subtract from the resisual at node 1.

        ;;      res_u(j) = res_u(j) - f(iu)
        ;;      res_p(j) = res_p(j) - f(ip)

        (let [j 0
              uL 0.0
              pL (nc/entry p j)
              uR (nc/entry u j)
              pR (nc/entry p j)
              fiu   (* 1/2 (+ (* nu pR) (* nu pL)
                              (/ nu (* Lr (- uR uL)))))
              fip (* 1/2 (+ (/ (+ uR    uL) Tr)
                            (/ nu (* Lr (- pR pL)))))

              new_res_u_j  (- (nc/entry res_u j) fiu)
              new_res_p_j  (- (nc/entry res_u j) fip)
              ]
          (do
            (clojure.pprint/pprint u)
            (clojure.pprint/pprint p)
            (clojure.pprint/pprint f)
            (clojure.pprint/pprint res_u)
            (clojure.pprint/pprint res_p)
            )
          (nc/entry! f iu fiu)
          (nc/entry! f ip fip)
          (nc/entry! res_u j new_res_u_j)
          (nc/entry! res_p j new_res_p_j))
        ;;   !-------------------------------------------
        ;;   ! Right boundary
        ;;   !-------------------------------------------

        ;;       j = nnodes

        ;;      uL = u(j)
        ;;      pL = p(j)

        ;;      uR = 0.0  !<- Boundary condition: u(1)=0.
        ;;      pR = pL   !<- No boundary condition: copy from the left.

        ;;     !Compute the numerical flux (upwind flux)

        ;;       f(iu) = half*( ( nu*pR + nu*pL )    + nu/Lr*( uR - uL ) )
        ;;       f(ip) = half*( (    uR +    uL )/Tr + nu/Lr*( pR - pL ) )

        ;;     !Add it to the resisual at node nnodes.

        ;;      res_u(j) = res_u(j) + f(iu)
        ;;      res_p(j) = res_p(j) + f(ip)

        (let [j (dec nnodes)
              uL (nc/entry u j)
              pL (nc/entry p j)
              uR 0.0
              pR uL
              fiu (* 1/2 (+ (* nu pR) (* nu pL)
                            (/ nu (* Lr (- uR uL)))))
              fip (* 1/2 (+ (/ (+ uR    uL) Tr)
                            (/ nu (* Lr (- pR pL)))))

              new_res_u_j  (+ (nc/entry res_u j) fiu)
              new_res_p_j  (+ (nc/entry res_u j) fip)
              ]
          (do
            (clojure.pprint/pprint u)
            (clojure.pprint/pprint p)
            (clojure.pprint/pprint f)
            (clojure.pprint/pprint res_u)
            (clojure.pprint/pprint res_p)
            )
          (nc/entry! f iu fiu)
          (nc/entry! f ip fip)
          (nc/entry! res_u j new_res_u_j)
          (nc/entry! res_p j new_res_p_j))
        ;; !------------------------------------------------------------
        ;; ! Add source terms, and finish the residual calculation.
        ;; !------------------------------------------------------------
        (doseq [j (range nnodes)]
          (let [new_res_u (+ (nc/entry res_u j)
                             (* nu nm/pi nm/pi
                                (nm/sin (* nm/pi (nc/entry x j)) )
                                (nc/entry vols j)))
                new_res_p (- (nc/entry res_p j)
                             (/ (nc/entry p j) Tr
                                (nc/entry vols j)))]
            (println (nc/entry res_u))
            (println new_res_p)
            (clojure.pprint/pprint (nc/entry x j))

            (nc/entry! res_u new_res_u)
            (nc/entry! res_p new_res_p)))

        ;;     do j = 1, nnodes
        ;;      res_u(j) = res_u(j) + nu*pi*pi*sin(pi*x(j))*vol(j)
        ;;      res_p(j) = res_p(j) - p(j)/Tr*vol(j)
        ;;     end do



        ;;  !-------------------------------------------------------------------------------
        ;;  !-------------------------------------------------------------------------------
        ;;  !  Check convergence: see if res_u and res_p are small enough.
        ;;  !-------------------------------------------------------------------------------
        ;;  !-------------------------------------------------------------------------------

        ;; !     Compute the maximum nodal residual (to check the convergence)
        ;;       res_max = max( maxval(abs(res_u)), maxval(abs(res_p)) )


        (do
         (let [res_max (max (nc/amax res_u) (nc/amax res_p))]


            ;; !     Stop if the tolerance (say, 1.0e-08) is reached.
            ;;       if ( res_max < 1.0e-08_dp ) exit pseudo_time_loop

            ;; !     Display the max nodal residual at every 100 iterations
            ;;       if (mod(k,500) == 0) then
            ;;         write(*,'(a5,i10,a20,es12.5)') \"Itr =\", k, \"   Max(nodal res) = \", res_max
            ;;       endif

            ;;  !-------------------------------------------------------------------------------
            ;;  !-------------------------------------------------------------------------------
            ;;  ! Update the solution by the forward Euler scheme.
            ;;  !-------------------------------------------------------------------------------
            ;;  !-------------------------------------------------------------------------------

            ;;             u = u + (dtau/h)*res_u
            ;;             p = p + (dtau/h)*res_p

            (if (and
                 (> 1.0e-08 res_max)
                 (< pesudotime 10000000))
              (let []
                (nc/axpby! 1.0 u (/ dtau heff) res_u)
                (nc/axpby! 1.0 p (/ dtau heff) res_p)
                (when (= 0 (mod pesudotime 500))
                  (do
                    (clojure.pprint/pprint u)
                    (clojure.pprint/pprint p)
                    (clojure.pprint/pprint f)
                    (clojure.pprint/pprint res_u)
                    (clojure.pprint/pprint res_p)
                    (println res_max)))
                (recur
                 (inc pesudotime)
                 x
                 u
                 p
                 f
                 res_u
                 res_p))

              (do
                (clojure.pprint/pprint u)
                (clojure.pprint/pprint p)
                (clojure.pprint/pprint f)
                (clojure.pprint/pprint res_u)
                (clojure.pprint/pprint res_p)
                (println res_max))
              )))))))

(let [x (irregular-grid 10)
      vols (dual-volumes x)
      spacing (mesh-spacing vols 10)
      gradients (gradients x)]

  (upwind-baby
   x vols spacing gradients))


(count (dual-volumes (irregular-grid 10)))



(dual-volumes (irregular-grid 10))
