(ns fluids.gen
  (:require [clojure.test.check :as stc]
            [clojure.spec.alpha :as s]
            [clojure.spec.gen.alpha :as gen]
            [uncomplicate.neanderthal.core :as nc]
            [uncomplicate.neanderthal.native :as nn]
            [uncomplicate.fluokitten.core :as fc]))


(defn ordered-dv-gen [& {:keys [count max] :or {count 4 max 1.0}}]
  (gen/fmap
   nn/dv
   (gen/fmap
    #(reductions + %)
    (s/gen
     (s/coll-of
      (s/double-in :min 0.1 :max max :infinite? false :NaN? false)
      :count count
      :distinct true)))))
