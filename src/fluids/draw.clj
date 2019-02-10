(ns fluids.draw
  (:require
   [fluids.core :as core]
   [uncomplicate.neanderthal.core :as nc]
   [clojure2d.core :as cc]
   [clojure2d.color :as color]))

(def vol-to-color
  (:orange-blue color/gradient-presets) )

(defn draw-line
  [cvs scale start end y vol count]
  (let [scale-start (* scale start )
        scale-end (* scale end )]
    (-> cvs
        (cc/set-color (vol-to-color  (* count vol)))
        (cc/line scale-start y scale-end y )
        (cc/set-color 255 255 255)
        (cc/ellipse scale-start y 5 5)
        (cc/ellipse scale-end y 5 5))
    cvs))

(defn draw-lines
  [canvas wavy-line vols count]
   (reduce
    (fn [canvas [[start end] vol]]
      (draw-line canvas 200 start end 100 vol count)
      canvas)
    canvas
    (map vector
         (partition 2 1 wavy-line)
         vols)) )


(def global-state (core/irregular-grid 10))

(def global-vols (core/dual-volumes global-state))

(defn drawer
  [canvas window ^long frameno state]
  (-> canvas
      (cc/set-background 0 0 0)
      (cc/set-color 255 255 255)
      (draw-lines
       global-state
       global-vols
       (nc/dim global-state)
       )))


(def window (cc/show-window
             {:canvas (cc/canvas 200 200 :mid) ; create canvas with mid quality
              :window-name "ellipse"           ; name window
              :w 400                           ; size of window (twice as canvas)
              :h 400
              :fps 10
              :hint :mid ;; hint for drawing canvas on window, mid quality (affects scalling 200 -> 400)
              :draw-fn  drawer})) ;; draw callback funtion



(comment
  (cc/close-window
   window))


