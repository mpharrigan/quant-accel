mkdir -p logs
accelerator interact --set_beta=0.0
for i in `seq 1 30`
do
    if [ $i -eq 15 ]
    then
        accelerator interact --set_beta=2.0
    fi

    accelerator TMat &>logs/round-$i-para-0 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-1 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-2 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-3 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-4 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-5 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-6 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-7 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-8 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-9 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-10 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-11 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-12 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-13 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-14 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-15 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-16 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-17 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-18 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-19 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-20 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-21 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-22 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-23 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-24 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-25 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-26 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-27 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-28 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-29 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-30 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-31 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-32 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-33 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-34 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-35 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-36 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-37 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-38 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-39 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-40 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-41 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-42 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-43 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-44 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-45 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-46 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-47 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-48 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-49 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-50 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-51 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-52 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-53 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-54 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-55 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-56 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-57 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-58 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-59 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-60 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-61 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-62 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-63 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-64 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-65 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-66 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-67 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-68 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-69 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-70 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-71 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-72 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-73 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-74 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-75 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-76 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-77 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-78 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-79 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-80 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-81 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-82 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-83 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-84 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-85 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-86 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-87 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-88 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-89 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-90 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-91 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-92 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-93 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-94 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-95 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-96 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-97 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-98 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-99 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-100 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-101 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-102 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-103 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-104 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-105 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-106 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-107 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-108 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-109 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-110 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-111 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-112 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-113 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-114 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-115 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-116 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-117 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-118 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-119 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-120 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-121 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-122 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-123 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-124 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-125 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-126 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-127 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-128 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-129 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-130 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-131 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-132 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-133 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-134 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-135 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-136 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-137 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-138 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-139 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-140 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-141 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-142 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-143 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-144 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-145 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-146 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-147 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-148 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-149 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-150 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-151 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-152 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-153 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-154 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-155 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-156 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-157 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-158 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-159 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-160 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-161 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-162 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-163 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-164 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-165 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-166 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-167 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-168 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-169 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-170 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-171 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-172 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-173 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-174 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-175 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-176 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-177 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-178 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-179 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-180 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-181 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-182 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-183 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-184 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-185 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-186 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-187 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-188 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-189 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-190 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-191 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-192 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-193 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-194 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-195 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-196 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-197 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-198 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-199 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-200 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-201 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-202 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-203 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-204 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-205 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-206 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-207 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-208 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-209 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-210 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-211 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-212 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-213 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-214 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-215 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-216 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-217 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-218 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-219 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-220 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-221 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-222 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-223 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-224 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-225 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-226 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-227 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-228 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-229 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-230 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-231 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-232 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-233 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-234 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-235 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-236 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-237 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-238 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-239 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-240 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-241 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-242 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-243 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-244 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-245 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-246 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-247 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-248 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-249 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-250 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-251 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-252 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-253 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-254 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-255 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-256 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-257 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-258 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-259 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-260 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-261 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-262 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-263 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-264 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-265 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-266 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-267 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-268 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-269 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-270 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-271 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-272 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-273 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-274 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-275 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-276 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-277 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-278 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-279 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-280 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-281 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-282 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-283 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-284 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-285 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-286 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-287 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-288 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-289 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-290 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-291 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-292 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-293 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-294 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-295 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-296 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-297 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-298 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-299 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-300 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-301 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-302 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-303 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-304 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-305 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-306 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-307 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-308 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-309 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-310 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-311 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-312 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-313 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-314 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-315 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-316 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-317 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-318 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-319 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-320 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-321 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-322 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-323 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-324 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-325 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-326 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-327 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-328 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-329 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-330 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-331 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-332 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-333 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-334 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-335 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-336 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-337 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-338 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-339 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-340 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-341 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-342 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-343 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-344 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-345 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-346 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-347 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-348 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-349 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-350 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-351 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-352 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-353 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-354 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-355 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-356 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-357 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-358 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-359 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-360 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-361 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-362 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-363 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-364 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-365 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-366 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-367 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-368 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-369 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-370 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-371 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-372 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-373 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-374 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-375 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-376 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-377 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-378 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-379 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-380 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-381 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-382 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-383 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-384 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-385 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-386 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-387 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-388 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-389 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-390 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-391 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-392 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-393 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-394 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-395 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-396 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-397 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-398 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-399 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-400 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-401 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-402 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-403 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-404 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-405 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-406 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-407 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-408 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-409 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-410 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-411 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-412 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-413 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-414 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-415 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-416 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-417 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-418 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-419 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-420 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-421 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-422 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-423 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-424 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-425 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-426 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-427 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-428 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-429 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-430 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-431 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-432 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-433 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-434 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-435 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-436 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-437 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-438 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-439 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-440 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-441 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-442 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-443 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-444 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-445 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-446 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-447 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-448 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-449 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-450 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-451 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-452 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-453 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-454 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-455 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-456 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-457 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-458 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-459 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-460 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-461 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-462 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-463 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-464 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-465 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-466 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-467 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-468 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-469 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-470 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-471 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-472 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-473 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-474 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-475 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-476 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-477 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-478 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-479 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-480 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-481 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-482 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-483 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-484 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-485 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-486 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-487 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-488 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-489 &
    sleep 10

    wait

    accelerator TMat &>logs/round-$i-para-490 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-491 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-492 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-493 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-494 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-495 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-496 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-497 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-498 &
    sleep 10

    accelerator TMat &>logs/round-$i-para-499 &
    sleep 10

    wait

    wait

    accelerator model &>logs/model-$i.log
done
